"""A transform that ignores `Nin` inputs and returns one constant value."""
struct ConstantTransform{Nin,T} <: AbstractTransform{Nin,1}
    value::T
end

ConstantTransform{Nin}(value) where {Nin} = Nin >= 0 ? ConstantTransform{Nin,typeof(value)}(value) : throw(ArgumentError("input dimension must be non-negative"))
ConstantTransform(value; n_inputs=0) = ConstantTransform{Int(n_inputs)}(value)

struct ConstantInverse{T<:ConstantTransform{0}} <: AbstractTransform{1,0}
    transform::T
end

@inline (transform::ConstantTransform{Nin,T})(::StaticVector{Nin}) where {Nin,T} = SVector{1,T}(transform.value)

function _comparison_tolerances(actual::Real, expected::Real, atol, rtol)
    T = promote_type(typeof(float(actual)), typeof(float(expected)))
    relative = rtol === nothing ? sqrt(eps(T)) : T(rtol)
    absolute = atol === nothing ? relative * max(one(T), abs(T(expected))) : T(atol)
    return absolute, relative
end

function _coordinates_match(actual::Real, expected::Real; atol=nothing, rtol=nothing)
    absolute, relative = _comparison_tolerances(actual, expected, atol, rtol)
    return isapprox(actual, expected; atol=absolute, rtol=relative)
end

inverse_transform(transform::ConstantTransform{0}) = ConstantInverse(transform)
inverse_transform(::ConstantTransform) = NoInverse()

@inline function (inverse::ConstantInverse)(coords::StaticVector{1})
    # A zero-dimensional result cannot carry NaN, so a mismatch must still throw.
    _coordinates_match(coords[1], inverse.transform.value) || throw(DomainError(coords, "coordinate does not match constant output $(inverse.transform.value)"))
    return SVector{0,_evaluation_type(eltype(coords))}()
end

inverse_transform(inverse::ConstantInverse) = inverse.transform

function transform_jacobian(transform::ConstantTransform{Nin}, coords::StaticVector{Nin}) where {Nin}
    T = Nin == 0 ? _evaluation_type(typeof(transform.value)) : promote_type(_evaluation_type(eltype(coords)), _evaluation_type(typeof(transform.value)))
    return SMatrix{1,Nin,T}(ntuple(_ -> zero(T), Nin))
end

function transform_jacobian(::ConstantInverse, coords::StaticVector{1})
    T = _evaluation_type(eltype(coords))
    return SMatrix{0,1,T}()
end

"""
    fix_inputs(transform, fixed; atol=nothing, rtol=nothing)

Return a lower-dimensional transform with selected inputs fixed to scalar
values. Programmatic integer indices are one-based; string and symbol keys
resolve through `input_names(transform)`.
"""
struct FixedInputsTransform{Nfree,Nfull,Nout,T<:AbstractTransform{Nfull,Nout},S,F,V,A,R} <: AbstractTransform{Nfree,Nout}
    transform::T
    sources::S
    free_indices::F
    fixed_values::V
    atol::A
    rtol::R
end

struct FixedInputsInverse{Nfree,Nfull,Nout,T<:FixedInputsTransform{Nfree,Nfull,Nout},I<:AbstractTransform{Nout,Nfull}} <: AbstractTransform{Nout,Nfree}
    transform::T
    inverse::I
end

function _validate_tolerance(value, name)
    value === nothing && return nothing
    value isa Real || throw(ArgumentError("$name must be a real number or nothing"))
    value >= 0 || throw(ArgumentError("$name must be non-negative"))
    return value
end

function _fixed_pairs(fixed)
    values = fixed isa Pair ? [fixed] : fixed isa NamedTuple || fixed isa AbstractDict ? collect(pairs(fixed)) : collect(fixed)
    isempty(values) && throw(ArgumentError("at least one input must be fixed"))
    all(value -> value isa Pair, values) || throw(ArgumentError("fixed inputs must be provided as pairs"))
    return values
end

function _fixed_index(key, names, count)
    if key isa Integer && !(key isa Bool)
        1 <= key <= count || throw(ArgumentError("fixed input index $key must be between 1 and $count"))
        return Int(key)
    elseif key isa AbstractString || key isa Symbol
        names === nothing && throw(ArgumentError("cannot fix input $key by name because the transform has no input names"))
        indices = findall(==(String(key)), names)
        isempty(indices) && throw(ArgumentError("fixed input name $key is not among $names"))
        length(indices) == 1 || throw(ArgumentError("fixed input name $key is ambiguous in $names"))
        return only(indices)
    end
    throw(ArgumentError("fixed input key $key must be an integer, string, or symbol"))
end

function fix_inputs(transform::AbstractTransform{Nfull,Nout}, fixed; atol=nothing, rtol=nothing) where {Nfull,Nout}
    input_pairs = _fixed_pairs(fixed)
    names = input_names(transform)
    indexed = map(input_pairs) do pair
        pair.second isa Real || throw(ArgumentError("fixed input values must be real scalars"))
        _fixed_index(pair.first, names, Nfull) => pair.second
    end
    indices = first.(indexed)
    length(unique(indices)) == length(indices) || throw(ArgumentError("the same input cannot be fixed more than once"))
    sort!(indexed; by=first)

    fixed_indices = Tuple(first.(indexed))
    raw_values = Tuple(last.(indexed))
    fixed_type = promote_type(map(value -> _evaluation_type(typeof(value)), raw_values)...)
    fixed_values = map(value -> fixed_type(value), raw_values)
    free_indices = Tuple(filter(index -> !(index in fixed_indices), 1:Nfull))
    sources = ntuple(Nfull) do index
        fixed_position = findfirst(==(index), fixed_indices)
        fixed_position === nothing ? findfirst(==(index), free_indices)::Int : -fixed_position
    end
    Nfree = length(free_indices)
    absolute = _validate_tolerance(atol, "atol")
    relative = _validate_tolerance(rtol, "rtol")
    return FixedInputsTransform{Nfree,Nfull,Nout,typeof(transform),typeof(sources),typeof(free_indices),typeof(fixed_values),typeof(absolute),typeof(relative)}(transform, sources, free_indices, fixed_values, absolute, relative)
end

function _fixed_evaluation_type(coords, fixed_values)
    coordinate_type = _evaluation_type(eltype(coords))
    fixed_types = map(value -> _evaluation_type(typeof(value)), fixed_values)
    return promote_type(coordinate_type, fixed_types...)
end

@inline function _full_inputs(transform::FixedInputsTransform{Nfree,Nfull}, coords::StaticVector{Nfree}) where {Nfree,Nfull}
    T = _fixed_evaluation_type(coords, transform.fixed_values)
    return SVector{Nfull,T}(ntuple(Nfull) do index
        source = transform.sources[index]
        source > 0 ? coords[source] : transform.fixed_values[-source]
    end)
end

@inline function (transform::FixedInputsTransform{Nfree})(coords::StaticVector{Nfree}) where {Nfree}
    return transform.transform(_full_inputs(transform, coords))
end

function inverse_transform(transform::FixedInputsTransform{Nfree,Nfull,Nout}) where {Nfree,Nfull,Nout}
    inverse = _inverse_transform_checked(transform.transform)
    inverse isa NoInverse && return inverse
    return FixedInputsInverse{Nfree,Nfull,Nout,typeof(transform),typeof(inverse)}(transform, inverse)
end

function _validate_fixed_inverse(inverse::FixedInputsInverse, coords, recovered)
    transform = inverse.transform
    for (index, source) in enumerate(transform.sources)
        source < 0 || continue
        expected = transform.fixed_values[-source]
        _coordinates_match(recovered[index], expected; atol=transform.atol, rtol=transform.rtol) || return map(value -> oftype(value, NaN), recovered)
    end
    roundtrip = transform.transform(recovered)
    # A coordinate off the fixed-input slice or outside the transform range has no preimage.
    all(index -> _coordinates_match(roundtrip[index], coords[index]; atol=transform.atol, rtol=transform.rtol), eachindex(coords)) || return map(value -> oftype(value, NaN), recovered)
    return recovered
end

@inline function _recover_fixed_inputs(inverse::FixedInputsInverse, coords)
    recovered = inverse.inverse(coords)
    return _validate_fixed_inverse(inverse, coords, recovered)
end

@inline function (inverse::FixedInputsInverse{Nfree})(coords::StaticVector) where {Nfree}
    recovered = _recover_fixed_inputs(inverse, coords)
    # With no free inputs the result cannot carry NaN, so an off-slice coordinate must throw.
    Nfree == 0 && any(isnan, recovered) && throw(DomainError(coords, "coordinate is outside the fixed-input slice"))
    T = _evaluation_type(eltype(recovered))
    return SVector{Nfree,T}(ntuple(index -> recovered[inverse.transform.free_indices[index]], Nfree))
end

inverse_transform(inverse::FixedInputsInverse) = inverse.transform

function compile_transform(transform::FixedInputsTransform)
    compiled = _compile_transform_checked(transform.transform)
    fixed_indices = filter(index -> transform.sources[index] < 0, eachindex(transform.sources))
    return fix_inputs(compiled, map(=>, fixed_indices, transform.fixed_values); atol=transform.atol, rtol=transform.rtol)
end

function transform_jacobian(transform::FixedInputsTransform{Nfree,Nfull,Nout}, coords::StaticVector{Nfree}) where {Nfree,Nfull,Nout}
    full = _full_inputs(transform, coords)
    jacobian = _transform_jacobian_checked(transform.transform, full)
    values = ntuple(Val(Nfree * Nout)) do index
        row = mod1(index, Nout)
        column = div(index - 1, Nout) + 1
        jacobian[row, transform.free_indices[column]]
    end
    return SMatrix{Nout,Nfree,eltype(jacobian)}(values)
end

function transform_jacobian(transform::FixedInputsTransform{0,Nfull,Nout}, coords::StaticVector{0}) where {Nfull,Nout}
    T = _fixed_evaluation_type(coords, transform.fixed_values)
    return SMatrix{Nout,0,T}()
end

function transform_jacobian(inverse::FixedInputsInverse{Nfree,Nfull,Nout}, coords::StaticVector{Nout}) where {Nfree,Nfull,Nout}
    _recover_fixed_inputs(inverse, coords)
    jacobian = _transform_jacobian_checked(inverse.inverse, coords)
    values = ntuple(Val(Nout * Nfree)) do index
        row = mod1(index, Nfree)
        column = div(index - 1, Nfree) + 1
        jacobian[inverse.transform.free_indices[row], column]
    end
    return SMatrix{Nfree,Nout,eltype(jacobian)}(values)
end

function transform_jacobian(inverse::FixedInputsInverse{0,Nfull,Nout}, coords::StaticVector{Nout}) where {Nfull,Nout}
    _recover_fixed_inputs(inverse, coords)
    return SMatrix{0,Nout,_evaluation_type(eltype(coords))}()
end

transform_name(transform::FixedInputsTransform) = transform_name(transform.transform)
output_names(transform::FixedInputsTransform) = output_names(transform.transform)

function input_names(transform::FixedInputsTransform{Nfree}) where {Nfree}
    names = input_names(transform.transform)
    names === nothing && return nothing
    return ntuple(index -> names[transform.free_indices[index]], Nfree)
end

transform_name(inverse::FixedInputsInverse) = transform_name(inverse.transform)
input_names(inverse::FixedInputsInverse) = output_names(inverse.transform)
output_names(inverse::FixedInputsInverse) = input_names(inverse.transform)
