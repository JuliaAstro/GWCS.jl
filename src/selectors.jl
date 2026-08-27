"""Non-coordinate mapping from numeric inputs to a region label."""
abstract type AbstractLabelMapper{N,L} end

function _selector_entries(values)
    entries = values isa AbstractDict ? collect(pairs(values)) : collect(values)
    all(entry -> entry isa Pair, entries) || throw(ArgumentError("selector entries must be pairs"))
    return entries
end

@generated function _apply_transform_at(transforms::T, index::Int, coords) where {T<:Tuple}
    expression = :(throw(BoundsError(transforms, index)))
    for branch in reverse(1:fieldcount(T))
        expression = :(index == $branch ? getfield(transforms, $branch)(coords) : $expression)
    end
    return expression
end

@generated function _transform_jacobian_at(transforms::T, index::Int, coords) where {T<:Tuple}
    expression = :(throw(BoundsError(transforms, index)))
    for branch in reverse(1:fieldcount(T))
        expression = :(index == $branch ? _transform_jacobian_checked(getfield(transforms, $branch), coords) : $expression)
    end
    return expression
end

input_dimension(::AbstractLabelMapper{N}) where {N} = N
label_type(::AbstractLabelMapper{N,L}) where {N,L} = L
input_names(::AbstractLabelMapper) = nothing
label_mapper_name(::AbstractLabelMapper) = nothing

function _label_dependencies(dependencies, ::Val{N}, context) where {N}
    dependencies === nothing && return SVector{N,Bool}(ntuple(_ -> true, N))
    dependencies isa AbstractVector || dependencies isa Tuple || throw(ArgumentError("$context must be a vector or tuple"))
    length(dependencies) == N || throw(DimensionMismatch("$context has length $(length(dependencies)), expected $N"))
    try
        return SVector{N,Bool}(dependencies)
    catch error
        error isa InexactError || error isa MethodError || rethrow()
        throw(ArgumentError("$context entries must be boolean-convertible"))
    end
end

label_dependencies(::AbstractLabelMapper{N}) where {N} = SVector{N,Bool}(ntuple(_ -> true, N))

function _label_dependencies_checked(mapper::AbstractLabelMapper{N}) where {N}
    return _label_dependencies(label_dependencies(mapper), Val(N), "label_dependencies(::$(typeof(mapper)))")
end

compile_label_mapper(mapper::AbstractLabelMapper) = mapper

function _compile_label_mapper_checked(mapper::AbstractLabelMapper)
    compiled = compile_label_mapper(mapper)
    compiled isa AbstractLabelMapper || throw(ArgumentError("compile_label_mapper(::$(typeof(mapper))) returned $(typeof(compiled)), expected AbstractLabelMapper"))
    input_dimension(compiled) == input_dimension(mapper) || throw(DimensionMismatch("compiled label mapper has input dimension $(input_dimension(compiled)), expected $(input_dimension(mapper))"))
    label_type(compiled) == label_type(mapper) || throw(ArgumentError("compiled label mapper has label type $(label_type(compiled)), expected $(label_type(mapper))"))
    return compiled
end

inverse_label_mapper(::AbstractLabelMapper) = NoInverse()

function _inverse_label_mapper_checked(mapper::AbstractLabelMapper)
    inverse = inverse_label_mapper(mapper)
    inverse isa NoInverse && return inverse
    inverse isa AbstractLabelMapper || throw(ArgumentError("inverse_label_mapper(::$(typeof(mapper))) returned $(typeof(inverse)), expected AbstractLabelMapper or NoInverse"))
    label_type(inverse) == label_type(mapper) || throw(ArgumentError("inverse label mapper has label type $(label_type(inverse)), expected $(label_type(mapper))"))
    return inverse
end

@inline map_label(mapper::AbstractLabelMapper{N}, coords::StaticVector{N}) where {N} = mapper(coords)

function map_label(mapper::AbstractLabelMapper{N}, coords::AbstractVector) where {N}
    length(coords) == N || throw(DimensionMismatch("label-mapper coordinate has length $(length(coords)), expected $N"))
    T = _evaluation_type(eltype(coords))
    return map_label(mapper, SVector{N,T}(coords))
end

function map_label(mapper::AbstractLabelMapper{N}, coords::Tuple{Vararg{Real,N}}) where {N}
    T = isempty(coords) ? Float64 : _evaluation_type(promote_type(map(typeof, coords)...))
    return map_label(mapper, SVector{N,T}(coords))
end

map_label(mapper::AbstractLabelMapper, coords::Real...) = map_label(mapper, coords)
map_label(mapper::AbstractLabelMapper{0}) = map_label(mapper, SVector{0,Float64}())

function _mapper_input_names(inputs, count)
    return _axis_names(inputs, count, "label-mapper input")
end

function _input_mapping(mapping, ::Val{Nin}, ::Val{Nout}) where {Nin,Nout}
    values = mapping === nothing ? ntuple(identity, Nout) : Tuple(Int.(mapping))
    length(values) == Nout || throw(DimensionMismatch("input mapping has length $(length(values)), expected $Nout"))
    all(index -> 1 <= index <= Nin, values) || throw(ArgumentError("input mapping indices must be between 1 and $Nin"))
    length(unique(values)) == Nout || throw(ArgumentError("input mapping indices must be unique"))
    return values
end

@inline function _project_inputs(coords::StaticVector, mapping::NTuple{N,Int}) where {N}
    return SVector{N}(ntuple(index -> coords[mapping[index]], N))
end

"""Attach a compatible inverse label mapper."""
struct WithInverseLabelMapper{Nin,Nout,L,F<:AbstractLabelMapper{Nin,L},I<:AbstractLabelMapper{Nout,L}} <: AbstractLabelMapper{Nin,L}
    forward::F
    backward::I
end

@inline (mapper::WithInverseLabelMapper)(coords::StaticVector) = map_label(mapper.forward, coords)
inverse_label_mapper(mapper::WithInverseLabelMapper) = WithInverseLabelMapper(mapper.backward, mapper.forward)
compile_label_mapper(mapper::WithInverseLabelMapper) = WithInverseLabelMapper(_compile_label_mapper_checked(mapper.forward), _compile_label_mapper_checked(mapper.backward))
label_dependencies(mapper::WithInverseLabelMapper) = _label_dependencies_checked(mapper.forward)
input_names(mapper::WithInverseLabelMapper) = input_names(mapper.forward)
label_mapper_name(mapper::WithInverseLabelMapper) = label_mapper_name(mapper.forward)

"""Two-dimensional array-backed region labels."""
struct ArrayLabelMapper{N,L,A<:AbstractMatrix{L},M<:NTuple{2,Int},NL,IN,NM} <: AbstractLabelMapper{N,L}
    labels::A
    mapping::M
    no_label::NL
    inputs::IN
    name::NM
end

function ArrayLabelMapper{N}(labels::AbstractMatrix{L}; inputs_mapping=nothing, no_label=nothing, inputs=nothing, name=nothing) where {N,L}
    mapping = _input_mapping(inputs_mapping, Val(N), Val(2))
    input_values = _mapper_input_names(inputs, N)
    name_value = name === nothing ? nothing : String(name)
    name_value === nothing || !isempty(name_value) || throw(ArgumentError("label-mapper name must not be empty"))
    return ArrayLabelMapper{N,L,typeof(labels),typeof(mapping),typeof(no_label),typeof(input_values),typeof(name_value)}(labels, mapping, no_label, input_values, name_value)
end

ArrayLabelMapper(labels::AbstractMatrix; kwargs...) = ArrayLabelMapper{2}(labels; kwargs...)

@inline function (mapper::ArrayLabelMapper{N})(coords::StaticVector{N}) where {N}
    projected = _project_inputs(coords, mapper.mapping)
    all(isfinite, projected) || return nothing
    rounded_x = floor(projected[1] + 0.5)
    rounded_y = floor(projected[2] + 0.5)
    0 <= rounded_x < size(mapper.labels, 2) || return nothing
    0 <= rounded_y < size(mapper.labels, 1) || return nothing
    x = firstindex(mapper.labels, 2) + Int(rounded_x)
    y = firstindex(mapper.labels, 1) + Int(rounded_y)
    label = @inbounds mapper.labels[y, x]
    return mapper.no_label !== nothing && isequal(label, mapper.no_label) ? nothing : label
end

function label_dependencies(mapper::ArrayLabelMapper{N}) where {N}
    return SVector{N,Bool}(ntuple(index -> index in mapper.mapping, N))
end

input_names(mapper::ArrayLabelMapper) = mapper.inputs
label_mapper_name(mapper::ArrayLabelMapper) = mapper.name

"""Callable region-label mapping with declared dimensions and dependencies."""
struct CallableLabelMapper{N,L,F,D,IN,NM} <: AbstractLabelMapper{N,L}
    callable::F
    dependencies::D
    inputs::IN
    name::NM
end

function CallableLabelMapper{N,L}(callable; dependencies=nothing, inputs=nothing, name=nothing) where {N,L}
    dependency_values = _label_dependencies(dependencies, Val(N), "callable label dependencies")
    input_values = _mapper_input_names(inputs, N)
    name_value = name === nothing ? nothing : String(name)
    name_value === nothing || !isempty(name_value) || throw(ArgumentError("label-mapper name must not be empty"))
    return CallableLabelMapper{N,L,typeof(callable),typeof(dependency_values),typeof(input_values),typeof(name_value)}(callable, dependency_values, input_values, name_value)
end

CallableLabelMapper(callable, nin::Integer, label_type::Type; kwargs...) = CallableLabelMapper{Int(nin),label_type}(callable; kwargs...)

@inline function (mapper::CallableLabelMapper{N,L})(coords::StaticVector{N}) where {N,L}
    result = mapper.callable(coords...)
    result === nothing && return nothing
    return convert(L, result)
end

label_dependencies(mapper::CallableLabelMapper) = mapper.dependencies
input_names(mapper::CallableLabelMapper) = mapper.inputs
label_mapper_name(mapper::CallableLabelMapper) = mapper.name

"""Adapt a one-output numeric transform as a label mapper."""
struct TransformLabelMapper{N,L,T<:AbstractTransform{N,1},NL,IN,NM} <: AbstractLabelMapper{N,L}
    transform::T
    no_label::NL
    inputs::IN
    name::NM
end

function TransformLabelMapper{L}(transform::AbstractTransform{N,1}; no_label=nothing, inputs=input_names(transform), name=transform_name(transform)) where {N,L}
    input_values = _mapper_input_names(inputs, N)
    name_value = name === nothing ? nothing : String(name)
    return TransformLabelMapper{N,L,typeof(transform),typeof(no_label),typeof(input_values),typeof(name_value)}(transform, no_label, input_values, name_value)
end

TransformLabelMapper(transform::AbstractTransform{N,1}; kwargs...) where {N} = TransformLabelMapper{Float64}(transform; kwargs...)

@inline function (mapper::TransformLabelMapper{N,L})(coords::StaticVector{N}) where {N,L}
    value = mapper.transform(coords)[1]
    mapper.no_label !== nothing && isequal(value, mapper.no_label) && return nothing
    return convert(L, value)
end

compile_label_mapper(mapper::TransformLabelMapper{N,L}) where {N,L} = TransformLabelMapper{L}(_compile_transform_checked(mapper.transform); no_label=mapper.no_label, inputs=mapper.inputs, name=mapper.name)
label_dependencies(mapper::TransformLabelMapper) = SVector(vec(_transform_dependencies_checked(mapper.transform)))
input_names(mapper::TransformLabelMapper) = mapper.inputs
label_mapper_name(mapper::TransformLabelMapper) = mapper.name

function _mapper_transforms(entries, ::Val{N}) where {N}
    values = Tuple(entry.second for entry in entries)
    all(transform -> transform isa AbstractTransform{N,1}, values) || throw(DimensionMismatch("all label-mapper transforms must have dimensions $N -> 1"))
    return values
end

function _mapper_label_type(::Type{L}) where {L}
    isconcretetype(L) || throw(ArgumentError("label type must be concrete"))
    return L
end

"""Tolerance-based numeric-key label mapping."""
struct DictionaryLabelMapper{N,L,K,T<:Tuple,M<:NTuple{1,Int},A,R,IN,NM} <: AbstractLabelMapper{N,L}
    keys::Vector{K}
    transforms::T
    mapping::M
    atol::A
    rtol::R
    inputs::IN
    name::NM
end

function DictionaryLabelMapper{N,L}(mapper; inputs_mapping=nothing, atol=1e-8, rtol=1e-5, inputs=nothing, name=nothing) where {N,L}
    _mapper_label_type(L)
    entries = _selector_entries(mapper)
    isempty(entries) && throw(ArgumentError("a dictionary label mapper requires at least one entry"))
    all(entry -> entry.first isa Real && isfinite(entry.first), entries) || throw(ArgumentError("dictionary label-mapper keys must be finite real numbers"))
    K = promote_type(map(entry -> typeof(entry.first), entries)...)
    keys = K[entry.first for entry in entries]
    transforms = _mapper_transforms(entries, Val(N))
    mapping = _input_mapping(inputs_mapping, Val(N), Val(1))
    atol >= 0 || throw(ArgumentError("atol must be nonnegative"))
    rtol >= 0 || throw(ArgumentError("rtol must be nonnegative"))
    input_values = _mapper_input_names(inputs, N)
    name_value = name === nothing ? nothing : String(name)
    return DictionaryLabelMapper{N,L,K,typeof(transforms),typeof(mapping),typeof(atol),typeof(rtol),typeof(input_values),typeof(name_value)}(keys, transforms, mapping, atol, rtol, input_values, name_value)
end

DictionaryLabelMapper(mapper, nin::Integer; label_type=Float64, kwargs...) = DictionaryLabelMapper{Int(nin),label_type}(mapper; kwargs...)

@inline function _matching_key(keys, value, atol, rtol)
    result = 0
    @inbounds for index in eachindex(keys)
        # This is the numpy isclose() comparison and is used for equivalence 
        # with the python implementation; isapprox would be stylistically better but is not equivalent.
        abs(value - keys[index]) <= atol + rtol * abs(value) && (result = index)
    end
    return result
end

@inline function (mapper::DictionaryLabelMapper{N,L})(coords::StaticVector{N}) where {N,L}
    key = coords[only(mapper.mapping)]
    isfinite(key) || return nothing
    index = _matching_key(mapper.keys, key, mapper.atol, mapper.rtol)
    iszero(index) && return nothing
    return convert(L, _apply_transform_at(mapper.transforms, index, coords)[1])
end

function compile_label_mapper(mapper::DictionaryLabelMapper{N,L}) where {N,L}
    transforms = map(_compile_transform_checked, mapper.transforms)
    entries = ntuple(index -> mapper.keys[index] => transforms[index], length(transforms))
    return DictionaryLabelMapper{N,L}(entries; inputs_mapping=mapper.mapping, atol=mapper.atol, rtol=mapper.rtol, inputs=mapper.inputs, name=mapper.name)
end

function label_dependencies(mapper::DictionaryLabelMapper{N}) where {N}
    result = MVector{N,Bool}(ntuple(index -> index in mapper.mapping, N))
    for transform in mapper.transforms
        result .|= SVector(vec(_transform_dependencies_checked(transform)))
    end
    return SVector(result)
end

input_names(mapper::DictionaryLabelMapper) = mapper.inputs
label_mapper_name(mapper::DictionaryLabelMapper) = mapper.name

"""Open-interval numeric-range label mapping."""
struct RangeLabelMapper{N,L,K,T<:Tuple,M<:NTuple{1,Int},IN,NM} <: AbstractLabelMapper{N,L}
    intervals::Vector{Tuple{K,K}}
    transforms::T
    mapping::M
    inputs::IN
    name::NM
end

function RangeLabelMapper{N,L}(mapper; inputs_mapping=nothing, inputs=nothing, name=nothing) where {N,L}
    _mapper_label_type(L)
    entries = _selector_entries(mapper)
    isempty(entries) && throw(ArgumentError("a range label mapper requires at least one entry"))
    all(entry -> entry.first isa Tuple && length(entry.first) == 2 && all(value -> value isa Real && isfinite(value), entry.first), entries) || throw(ArgumentError("range label-mapper keys must be finite two-value tuples"))
    K = promote_type((typeof(value) for entry in entries for value in entry.first)...)
    intervals = Tuple{K,K}[(K(entry.first[1]), K(entry.first[2])) for entry in entries]
    all(interval -> interval[1] < interval[2], intervals) || throw(ArgumentError("range label-mapper lower limits must be less than upper limits"))
    permutation = sortperm(intervals; by=first)
    intervals = intervals[permutation]
    transforms_unsorted = _mapper_transforms(entries, Val(N))
    transforms = ntuple(index -> transforms_unsorted[permutation[index]], length(permutation))
    for index in 1:(length(intervals) - 1)
        intervals[index][2] > intervals[index + 1][1] && throw(ArgumentError("range label-mapper intervals must not overlap"))
    end
    mapping = _input_mapping(inputs_mapping, Val(N), Val(1))
    input_values = _mapper_input_names(inputs, N)
    name_value = name === nothing ? nothing : String(name)
    return RangeLabelMapper{N,L,K,typeof(transforms),typeof(mapping),typeof(input_values),typeof(name_value)}(intervals, transforms, mapping, input_values, name_value)
end

RangeLabelMapper(mapper, nin::Integer; label_type=Float64, kwargs...) = RangeLabelMapper{Int(nin),label_type}(mapper; kwargs...)

@inline function _matching_interval(intervals, value)
    @inbounds for index in eachindex(intervals)
        lower, upper = intervals[index]
        lower < value < upper && return index
    end
    return 0
end

@inline function (mapper::RangeLabelMapper{N,L})(coords::StaticVector{N}) where {N,L}
    key = coords[only(mapper.mapping)]
    isfinite(key) || return nothing
    index = _matching_interval(mapper.intervals, key)
    iszero(index) && return nothing
    return convert(L, _apply_transform_at(mapper.transforms, index, coords)[1])
end

function compile_label_mapper(mapper::RangeLabelMapper{N,L}) where {N,L}
    transforms = map(_compile_transform_checked, mapper.transforms)
    entries = ntuple(index -> mapper.intervals[index] => transforms[index], length(transforms))
    return RangeLabelMapper{N,L}(entries; inputs_mapping=mapper.mapping, inputs=mapper.inputs, name=mapper.name)
end

function label_dependencies(mapper::RangeLabelMapper{N}) where {N}
    result = MVector{N,Bool}(ntuple(index -> index in mapper.mapping, N))
    for transform in mapper.transforms
        result .|= SVector(vec(_transform_dependencies_checked(transform)))
    end
    return SVector(result)
end

input_names(mapper::RangeLabelMapper) = mapper.inputs
label_mapper_name(mapper::RangeLabelMapper) = mapper.name
