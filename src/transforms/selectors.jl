function _region_entries(selector, ::Val{Nin}, ::Val{Nout}) where {Nin,Nout}
    entries = _selector_entries(selector)
    isempty(entries) && throw(ArgumentError("a region selector requires at least one transform"))
    labels = map(entry -> entry.first, entries)
    label_types = unique(typeof.(labels))
    length(label_types) == 1 || throw(ArgumentError("region labels must have one concrete type"))
    L = only(label_types)
    isconcretetype(L) || throw(ArgumentError("region label type must be concrete"))
    transforms = Tuple(entry.second for entry in entries)
    all(transform -> transform isa AbstractTransform{Nin,Nout}, transforms) || throw(DimensionMismatch("all selected transforms must have dimensions $Nin -> $Nout"))
    lookup = Dict{L,Int}()
    for (index, label) in enumerate(labels)
        haskey(lookup, label) && throw(ArgumentError("region labels must be unique"))
        lookup[label] = index
    end
    return L, collect(L, labels), transforms, lookup
end

"""Discontinuous transform selected by a non-coordinate region label."""
struct RegionSelector{Nin,Nout,L,M<:AbstractLabelMapper{Nin,L},Ts<:Tuple,D<:AbstractDict{L,Int},U<:Real} <: AbstractTransform{Nin,Nout}
    label_mapper::M
    labels::Vector{L}
    transforms::Ts
    lookup::D
    undefined_value::U
end

function RegionSelector(label_mapper::AbstractLabelMapper{Nin,L}, selector; undefined_value=NaN) where {Nin,L}
    undefined_value isa Real || throw(ArgumentError("undefined_value must be a real number"))
    entries = _selector_entries(selector)
    isempty(entries) && throw(ArgumentError("a region selector requires at least one transform"))
    Nout = output_dimension(first(entries).second)
    label_type_value, labels, transforms, lookup = _region_entries(entries, Val(Nin), Val(Nout))
    label_type_value == L || throw(ArgumentError("region labels have type $label_type_value, expected $L from the label mapper"))
    return RegionSelector{Nin,Nout,L,typeof(label_mapper),typeof(transforms),typeof(lookup),typeof(undefined_value)}(label_mapper, labels, transforms, lookup, undefined_value)
end

@inline function _undefined_region_output(transform::RegionSelector{Nin,Nout}, coords::StaticVector{Nin}) where {Nin,Nout}
    T = _evaluation_type(eltype(coords))
    return SVector{Nout,T}(ntuple(_ -> T(transform.undefined_value), Nout))
end

@inline function (transform::RegionSelector{Nin,Nout})(coords::StaticVector{Nin}) where {Nin,Nout}
    label = map_label(transform.label_mapper, coords)
    label === nothing && return _undefined_region_output(transform, coords)
    index = get(transform.lookup, label, 0)
    iszero(index) && return _undefined_region_output(transform, coords)
    return _apply_transform_at(transform.transforms, index, coords)
end

function selected_transform(transform::RegionSelector, label)
    index = get(transform.lookup, label, 0)
    iszero(index) && throw(KeyError(label))
    return transform.transforms[index]
end

function compile_transform(transform::RegionSelector)
    mapper = _compile_label_mapper_checked(transform.label_mapper)
    transforms = map(_compile_transform_checked, transform.transforms)
    selector = ntuple(index -> transform.labels[index] => transforms[index], length(transforms))
    return RegionSelector(mapper, selector; undefined_value=transform.undefined_value)
end

function inverse_transform(transform::RegionSelector)
    mapper = _inverse_label_mapper_checked(transform.label_mapper)
    mapper isa NoInverse && return mapper
    inverses = map(_inverse_transform_checked, transform.transforms)
    any(inverse -> inverse isa NoInverse, inverses) && return NoInverse()
    input_dimension(mapper) == output_dimension(transform) || throw(DimensionMismatch("inverse label mapper has input dimension $(input_dimension(mapper)), expected $(output_dimension(transform))"))
    selector = ntuple(index -> transform.labels[index] => inverses[index], length(inverses))
    return RegionSelector(mapper, selector; undefined_value=transform.undefined_value)
end

function transform_jacobian(transform::RegionSelector{Nin,Nout}, coords::StaticVector{Nin}) where {Nin,Nout}
    # Outside every region the selector has no derivative.
    undefined = SMatrix{Nout,Nin}(ntuple(_ -> _evaluation_type(eltype(coords))(NaN), Nin * Nout))
    label = map_label(transform.label_mapper, coords)
    label === nothing && return undefined
    index = get(transform.lookup, label, 0)
    iszero(index) && return undefined
    return _transform_jacobian_at(transform.transforms, index, coords)
end

function transform_dependencies(transform::RegionSelector{Nin,Nout}) where {Nin,Nout}
    dependencies = MMatrix{Nout,Nin,Bool}(false for _ in 1:Nout, _ in 1:Nin)
    for branch in transform.transforms
        dependencies .|= _transform_dependencies_checked(branch)
    end
    selection = _label_dependencies_checked(transform.label_mapper)
    for row in 1:Nout, column in 1:Nin
        dependencies[row, column] |= selection[column]
    end
    return SMatrix(dependencies)
end
