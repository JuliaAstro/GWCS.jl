transform_dependencies(::IdentityTransform{N}) where {N} = SMatrix{N,N,Bool}(I)
transform_dependencies(transform::_ApproximateInverseTransform) = _transform_dependencies_checked(transform.forward)
transform_dependencies(::Shift) = SMatrix{1,1,Bool}(true)
transform_dependencies(transform::Scale) = SMatrix{1,1,Bool}(!iszero(transform.factor))

function transform_dependencies(transform::AffineTransform{Nin,Nout}) where {Nin,Nout}
    return SMatrix{Nout,Nin,Bool}(map(!iszero, transform.matrix))
end

function transform_dependencies(transform::RemapAxes{Nin,Nout}) where {Nin,Nout}
    values = ntuple(Val(Nin * Nout)) do index
        row = mod1(index, Nout)
        column = div(index - 1, Nout) + 1
        transform.mapping[row] == column
    end
    return SMatrix{Nout,Nin,Bool}(values)
end

function transform_dependencies(transform::Polynomial1D)
    depends = any(!iszero, Iterators.drop(transform.coefficients, 1))
    return SMatrix{1,1,Bool}(depends)
end

function transform_dependencies(transform::Polynomial2D{N}) where {N}
    depends_x = any(!iszero(transform.coefficients[i + 1, j + 1]) for j in 0:(N - 1) for i in 1:(N - 1 - j))
    depends_y = any(!iszero(transform.coefficients[i + 1, j + 1]) for j in 1:(N - 1) for i in 0:(N - 1 - j))
    return SMatrix{1,2,Bool}(depends_x, depends_y)
end

transform_dependencies(::SphericalToCartesian) = SMatrix{3,2,Bool}(true, true, false, true, true, true)
transform_dependencies(::CartesianToSpherical) = SMatrix{2,3,Bool}(true, true, true, true, false, true)
transform_dependencies(::IntermediateToNative) = _dense_dependencies(Val(2), Val(2))
transform_dependencies(::NativeToIntermediate) = _dense_dependencies(Val(2), Val(2))
transform_dependencies(::NativeToCelestial) = _dense_dependencies(Val(2), Val(2))
transform_dependencies(::CelestialToNative) = _dense_dependencies(Val(2), Val(2))

function transform_dependencies(transform::Rotation3D)
    return SMatrix{3,3,Bool}(map(!iszero, transform.matrix))
end

function transform_dependencies(::ConstantTransform{Nin}) where {Nin}
    return SMatrix{1,Nin,Bool}(ntuple(_ -> false, Nin))
end

transform_dependencies(::ConstantInverse) = SMatrix{0,1,Bool}()
transform_dependencies(transform::AnnotatedTransform) = _transform_dependencies_checked(transform.transform)
transform_dependencies(transform::WithInverse) = _transform_dependencies_checked(transform.forward)

function transform_dependencies(transform::FixedInputsTransform{Nfree,Nfull,Nout}) where {Nfree,Nfull,Nout}
    dependencies = _transform_dependencies_checked(transform.transform)
    values = ntuple(Val(Nfree * Nout)) do index
        row = mod1(index, Nout)
        column = div(index - 1, Nout) + 1
        dependencies[row, transform.free_indices[column]]
    end
    return SMatrix{Nout,Nfree,Bool}(values)
end

function transform_dependencies(transform::FixedInputsInverse{Nfree,Nfull,Nout}) where {Nfree,Nfull,Nout}
    dependencies = _transform_dependencies_checked(transform.inverse)
    values = ntuple(Val(Nout * Nfree)) do index
        row = mod1(index, Nfree)
        column = div(index - 1, Nfree) + 1
        dependencies[transform.transform.free_indices[row], column]
    end
    return SMatrix{Nfree,Nout,Bool}(values)
end

function _compose_dependencies(first::StaticMatrix{Nmid,Nin,Bool}, second::StaticMatrix{Nout,Nmid,Bool}) where {Nin,Nmid,Nout}
    values = ntuple(Val(Nin * Nout)) do index
        row = mod1(index, Nout)
        column = div(index - 1, Nout) + 1
        any(intermediate -> second[row, intermediate] && first[intermediate, column], 1:Nmid)
    end
    return SMatrix{Nout,Nin,Bool}(values)
end

_composed_dependencies(transforms::Tuple{Any}) = _transform_dependencies_checked(first(transforms))

function _composed_dependencies(transforms::Tuple)
    first_dependencies = _transform_dependencies_checked(first(transforms))
    remaining_dependencies = _composed_dependencies(Base.tail(transforms))
    return _compose_dependencies(first_dependencies, remaining_dependencies)
end

transform_dependencies(transform::ComposedTransform) = _composed_dependencies(transform.transforms)

@inline _parallel_dependency_entry(::Tuple{}, row, column, input_offset, output_offset) = false

@inline function _parallel_dependency_entry(transforms::Tuple, row, column, input_offset, output_offset)
    transform = first(transforms)
    dependencies = _transform_dependencies_checked(transform)
    Nout, Nin = size(dependencies)
    if output_offset < row <= output_offset + Nout && input_offset < column <= input_offset + Nin
        return dependencies[row - output_offset, column - input_offset]
    end
    return _parallel_dependency_entry(Base.tail(transforms), row, column, input_offset + Nin, output_offset + Nout)
end

function transform_dependencies(transform::ParallelTransform{Nin,Nout}) where {Nin,Nout}
    values = ntuple(Val(Nin * Nout)) do index
        row = mod1(index, Nout)
        column = div(index - 1, Nout) + 1
        _parallel_dependency_entry(transform.transforms, row, column, 0, 0)
    end
    return SMatrix{Nout,Nin,Bool}(values)
end

function transform_dependencies(transform::ArithmeticTransform{Nin,Nout}) where {Nin,Nout}
    left = _transform_dependencies_checked(transform.left)
    right = _transform_dependencies_checked(transform.right)
    return SMatrix{Nout,Nin,Bool}(ntuple(index -> left[index] || right[index], Nin * Nout))
end

function transform_dependencies(transform::CallableTransform{Nin,Nout}) where {Nin,Nout}
    transform.dependencies === nothing && return _dense_dependencies(Val(Nout), Val(Nin))
    return transform.dependencies
end

function transform_dependencies(transform::InverseCallableTransform{Nin,Nout}) where {Nin,Nout}
    transform.transform.inverse_dependencies === nothing && return _dense_dependencies(Val(Nout), Val(Nin))
    return transform.transform.inverse_dependencies
end

function _dependency_matrix(entry, ::Val{Nout}, ::Val{Nin}) where {Nout,Nin}
    values = ntuple(index -> entry(mod1(index, Nout), div(index - 1, Nout) + 1)::Bool, Val(Nin * Nout))
    return SMatrix{Nout,Nin,Bool}(values)
end

transform_dependencies(table::FITSLookupTable{N,Nin}) where {N,Nin} = _dependency_matrix((_, column) -> column in table.axes, Val(1), Val(Nin))

# Dependencies of the per-axis lookup offsets on the first two axes, optionally including the identity term.
function _lookup_stage_dependencies(tables, ::Val{N}, include_identity) where {N}
    return _dependency_matrix(Val(N), Val(N)) do row, column
        include_identity && row == column || row <= 2 && tables[row] !== nothing && column in tables[row].axes
    end
end

function transform_dependencies(transform::FITSLookupDistortion{N}) where {N}
    detector = _lookup_stage_dependencies(transform.detector_to_image, Val(N), true)
    corrected = _transform_dependencies_checked(transform.sip) .| _lookup_stage_dependencies(transform.cpdis, Val(N), false)
    return _compose_dependencies(detector, corrected)
end

_tpd_axis_dependency(::Nothing, row, column) = row == column
_tpd_axis_dependency(model::_FITSTPDAxis, row, column) = column in model.axes || model.docorr && row == column

transform_dependencies(transform::_FITSTPDTransform{N}) where {N} = _dependency_matrix((row, column) -> _tpd_axis_dependency(transform.models[row], row, column), Val(N), Val(N))

function transform_dependencies(transform::_FITSPriorDistortion{N}) where {N}
    detector = _lookup_stage_dependencies(transform.detector_to_image, Val(N), true)
    sip = _transform_dependencies_checked(transform.sip)
    lookup = _compose_dependencies(detector, sip .| _lookup_stage_dependencies(transform.cpdis, Val(N), false))
    tpd = _compose_dependencies(_compose_dependencies(detector, sip), _transform_dependencies_checked(transform.tpd))
    return _dependency_matrix((row, column) -> transform.tpd.models[row] === nothing ? lookup[row, column] : tpd[row, column], Val(N), Val(N))
end

function _tabular_dependency(tables, row, column)
    for table in tables
        row in table.axes && return column in table.axes
    end
    return row == column
end

transform_dependencies(transform::Union{FITSTabularTransform{N},FITSTabularInverse{N}}) where {N} = _dependency_matrix((row, column) -> _tabular_dependency(transform.tables, row, column), Val(N), Val(N))
