"""Input-domain bounds for an `N`-dimensional transform."""
abstract type AbstractBoundingBox{N} end

"""Rectangular bounds on selected input axes."""
struct BoundingBox{N,T<:Real,K,A<:NTuple{K,Int}} <: AbstractBoundingBox{N}
    intervals::NTuple{K,Tuple{T,T}}
    axes::A
    order::Symbol
end

function _bounding_box_values(intervals)
    values = Tuple(begin
        length(interval) == 2 || throw(ArgumentError("bounding-box intervals must contain exactly two bounds"))
        lower, upper = interval
        lower isa Real && upper isa Real || throw(ArgumentError("bounding-box bounds must be real numbers"))
        lower <= upper || throw(ArgumentError("bounding-box lower limits must not exceed upper limits"))
        promote(lower, upper)
    end for interval in intervals)
    isempty(values) && return (), Float64
    T = promote_type(map(interval -> promote_type(typeof(interval[1]), typeof(interval[2])), values)...)
    return ntuple(index -> (T(values[index][1]), T(values[index][2])), length(values)), T
end

function BoundingBox{N}(intervals, axes; order=:F) where {N}
    order in (:F, :C) || throw(ArgumentError("bounding-box order must be :F or :C"))
    axis_values = Tuple(Int(axis) for axis in axes)
    interval_values, T = _bounding_box_values(intervals)
    length(interval_values) == length(axis_values) || throw(DimensionMismatch("bounding box has $(length(interval_values)) intervals but $(length(axis_values)) bounded axes"))
    all(axis -> 1 <= axis <= N, axis_values) || throw(ArgumentError("bounded axes must be between 1 and $N"))
    length(unique(axis_values)) == length(axis_values) || throw(ArgumentError("bounded axes must be unique"))
    if order === :C
        interval_values = reverse(interval_values)
    end
    permutation = sortperm(collect(axis_values))
    canonical_axes = ntuple(index -> axis_values[permutation[index]], length(axis_values))
    canonical_intervals = ntuple(index -> interval_values[permutation[index]], length(interval_values))
    K = length(canonical_axes)
    return BoundingBox{N,T,K,typeof(canonical_axes)}(canonical_intervals, canonical_axes, :F)
end

BoundingBox(intervals::Tuple; order=:F) = BoundingBox{length(intervals)}(intervals, ntuple(identity, length(intervals)); order)
BoundingBox(intervals::AbstractVector; order=:F) = BoundingBox(Tuple(intervals); order)
BoundingBox{N,T}(intervals::NTuple{N,Tuple{T,T}}, order::Symbol) where {N,T<:Real} = BoundingBox(intervals; order)

input_dimension(::AbstractBoundingBox{N}) where {N} = N
bounded_axes(box::BoundingBox) = box.axes
ignored_axes(box::BoundingBox{N}) where {N} = Tuple(filter(axis -> !(axis in box.axes), 1:N))

function bounds(box::BoundingBox, axis::Integer)
    1 <= axis <= input_dimension(box) || throw(BoundsError(box, axis))
    index = findfirst(==(axis), box.axes)
    return index === nothing ? nothing : box.intervals[index]
end

@inline function Base.in(coords::StaticVector{N}, box::BoundingBox{N}) where {N}
    return all(index -> begin
        lower, upper = box.intervals[index]
        lower <= coords[box.axes[index]] <= upper
    end, eachindex(box.axes))
end

function Base.in(coords::AbstractVector, box::AbstractBoundingBox{N}) where {N}
    length(coords) == N || throw(DimensionMismatch("coordinate has length $(length(coords)), expected $N"))
    return SVector{N}(coords) in box
end

struct _BoundingBoxLookup{K,V} <: AbstractDict{K,V}
    entries::Vector{Pair{K,V}}
end

Base.length(lookup::_BoundingBoxLookup) = length(lookup.entries)
Base.iterate(lookup::_BoundingBoxLookup, state...) = iterate(lookup.entries, state...)

function _lookup_bounding_box_index(lookup::_BoundingBoxLookup, key)
    first_index = firstindex(lookup.entries)
    last_index = lastindex(lookup.entries)
    while first_index <= last_index
        middle = (first_index + last_index) >>> 1
        entry = @inbounds lookup.entries[middle]
        isequal(entry.first, key) && return middle
        if isless(entry.first, key)
            first_index = middle + 1
        else
            last_index = middle - 1
        end
    end
    return 0
end

function Base.getindex(lookup::_BoundingBoxLookup, key)
    index = _lookup_bounding_box_index(lookup, key)
    iszero(index) && throw(KeyError(key))
    return @inbounds lookup.entries[index].second
end

function Base.get(lookup::_BoundingBoxLookup, key, default)
    index = _lookup_bounding_box_index(lookup, key)
    return iszero(index) ? default : @inbounds(lookup.entries[index].second)
end
Base.haskey(lookup::_BoundingBoxLookup, key) = !iszero(_lookup_bounding_box_index(lookup, key))

"""Selector-keyed collection of rectangular input domains."""
struct CompoundBoundingBox{N,K,S<:NTuple{K,Int},G<:NTuple{K,Bool},I<:Tuple,D<:AbstractDict} <: AbstractBoundingBox{N}
    selector_axes::S
    selector_ignored::G
    globally_ignored::I
    boxes::D
end

function _compound_key(selector, count)
    values = selector isa Tuple ? selector : selector isa Real ? (selector,) : Tuple(selector)
    length(values) == count || throw(DimensionMismatch("selector has length $(length(values)), expected $count"))
    all(value -> value isa Real, values) || throw(ArgumentError("selector values must be real numbers"))
    return values
end

function _normalize_compound_boxes(boxes, ::Val{N}, ::Val{K}) where {N,K}
    entries = collect(pairs(boxes))
    isempty(entries) && throw(ArgumentError("a compound bounding box requires at least one entry"))
    all(entry -> entry.second isa BoundingBox{N}, entries) || throw(DimensionMismatch("all compound bounding boxes must have input dimension $N"))
    axes = bounded_axes(first(entries).second)
    all(entry -> bounded_axes(entry.second) == axes, entries) || throw(ArgumentError("compound bounding-box entries must constrain the same axes"))
    key_types = ntuple(index -> promote_type(Float64, map(entry -> typeof(_compound_key(entry.first, K)[index]), entries)...), K)
    Key = Tuple{key_types...}
    bound_types = map(entries) do entry
        box = entry.second
        isempty(box.intervals) ? Float64 : eltype(first(box.intervals))
    end
    T = promote_type(bound_types...)
    Box = BoundingBox{N}(ntuple(index -> begin
        lower, upper = first(entries).second.intervals[index]
        (T(lower), T(upper))
    end, length(axes)), axes)
    result = Pair{Key,typeof(Box)}[]
    seen = Set{Key}()
    for entry in entries
        key_values = _compound_key(entry.first, K)
        key = Key(ntuple(index -> key_types[index](key_values[index]), K))
        key in seen && throw(ArgumentError("compound bounding-box selector keys must be unique"))
        push!(seen, key)
        box = BoundingBox{N}(ntuple(index -> begin
            lower, upper = entry.second.intervals[index]
            (T(lower), T(upper))
        end, length(axes)), axes)
        push!(result, key => box)
    end
    sort!(result; by=first)
    return _BoundingBoxLookup(result)
end

function CompoundBoundingBox{N}(boxes, selector_axes; selector_ignored=true, ignored=()) where {N}
    selectors = Tuple(Int(axis) for axis in selector_axes)
    isempty(selectors) && throw(ArgumentError("a compound bounding box requires at least one selector axis"))
    all(axis -> 1 <= axis <= N, selectors) || throw(ArgumentError("selector axes must be between 1 and $N"))
    length(unique(selectors)) == length(selectors) || throw(ArgumentError("selector axes must be unique"))
    K = length(selectors)
    if selector_ignored isa Bool
        selector_flags = ntuple(_ -> selector_ignored, K)
    else
        all(value -> value isa Bool, selector_ignored) || throw(ArgumentError("selector ignore flags must be booleans"))
        selector_flags = Tuple(selector_ignored)
    end
    length(selector_flags) == K || throw(DimensionMismatch("selector ignore flags have length $(length(selector_flags)), expected $K"))
    global_axes = Tuple(Int(axis) for axis in ignored)
    all(axis -> 1 <= axis <= N, global_axes) || throw(ArgumentError("ignored axes must be between 1 and $N"))
    length(unique(global_axes)) == length(global_axes) || throw(ArgumentError("ignored axes must be unique"))
    isempty(intersect(selectors, global_axes)) || throw(ArgumentError("selector axes cannot also be globally ignored"))
    normalized = _normalize_compound_boxes(boxes, Val(N), Val(K))
    constrained = bounded_axes(first(values(normalized)))
    for (index, axis) in enumerate(selectors)
        selector_flags[index] == !(axis in constrained) || throw(ArgumentError("selector axis $axis ignore flag is inconsistent with the selected bounding boxes"))
    end
    isempty(intersect(constrained, global_axes)) || throw(ArgumentError("globally ignored axes cannot also be bounded"))
    accounted = union(Set(constrained), Set(selectors), Set(global_axes))
    accounted == Set(1:N) || throw(ArgumentError("every input axis must be bounded, selected, or globally ignored"))
    return CompoundBoundingBox{N,K,typeof(selectors),typeof(selector_flags),typeof(global_axes),typeof(normalized)}(selectors, selector_flags, global_axes, normalized)
end

function CompoundBoundingBox(boxes, selector_axes; kwargs...)
    isempty(boxes) && throw(ArgumentError("a compound bounding box requires at least one entry"))
    return CompoundBoundingBox{input_dimension(first(values(boxes)))}(boxes, selector_axes; kwargs...)
end

selector_axes(box::CompoundBoundingBox) = box.selector_axes
selector_ignored(box::CompoundBoundingBox) = box.selector_ignored
ignored_axes(box::CompoundBoundingBox) = box.globally_ignored

function _selected_bounding_box_index(box::CompoundBoundingBox{N,K}, selector) where {N,K}
    key = convert(keytype(box.boxes), _compound_key(selector, K))
    return _lookup_bounding_box_index(box.boxes, key)
end

function select_bounding_box(box::CompoundBoundingBox, selector)
    index = _selected_bounding_box_index(box, selector)
    iszero(index) && throw(KeyError(selector))
    return @inbounds box.boxes.entries[index].second
end

Base.getindex(box::CompoundBoundingBox, selector) = select_bounding_box(box, selector)

@inline function Base.in(coords::StaticVector{N}, box::CompoundBoundingBox{N,K}) where {N,K}
    selector = ntuple(index -> coords[box.selector_axes[index]], K)
    index = _selected_bounding_box_index(box, selector)
    iszero(index) && return false
    return coords in @inbounds(box.boxes.entries[index].second)
end

function _free_axis_map(count, fixed_axes)
    free = Tuple(filter(axis -> !(axis in fixed_axes), 1:count))
    return free, axis -> findfirst(==(axis), free)
end

function _reduce_bounding_box(box::BoundingBox{N}, fixed) where {N}
    fixed_axes = Tuple(first.(fixed))
    free, new_axis = _free_axis_map(N, fixed_axes)
    for (axis, value) in fixed
        interval = bounds(box, axis)
        interval === nothing && continue
        interval[1] <= value <= interval[2] || throw(DomainError(value, "fixed input $axis lies outside the WCS bounding box"))
    end
    retained = Tuple(index for index in eachindex(box.axes) if !(box.axes[index] in fixed_axes))
    isempty(retained) && return nothing
    axes = map(index -> new_axis(box.axes[index])::Int, retained)
    intervals = map(index -> box.intervals[index], retained)
    return BoundingBox{length(free)}(intervals, axes)
end

function _reduce_bounding_box(box::CompoundBoundingBox{N,K}, fixed) where {N,K}
    fixed_map = Dict(fixed)
    selector_fixed = ntuple(index -> haskey(fixed_map, box.selector_axes[index]), K)
    surviving = Pair[]
    for (key, selected) in box.boxes
        all(index -> !selector_fixed[index] || isequal(key[index], fixed_map[box.selector_axes[index]]), 1:K) || continue
        reduced = try
            _reduce_bounding_box(selected, fixed)
        catch error
            error isa DomainError || rethrow()
            continue
        end
        remaining_key = Tuple(key[index] for index in 1:K if !selector_fixed[index])
        push!(surviving, remaining_key => reduced)
    end
    isempty(surviving) && throw(DomainError(values(fixed_map), "fixed selector inputs do not match a compound bounding box"))
    remaining_indices = Tuple(index for index in 1:K if !selector_fixed[index])
    isempty(remaining_indices) && return last(only(surviving))
    surviving = map(surviving) do entry
        reduced = last(entry)
        first(entry) => (reduced === nothing ? BoundingBox{N - length(fixed)}((), ()) : reduced)
    end
    fixed_axes = Tuple(keys(fixed_map))
    free, new_axis = _free_axis_map(N, fixed_axes)
    selectors = map(index -> new_axis(box.selector_axes[index])::Int, remaining_indices)
    flags = map(index -> box.selector_ignored[index], remaining_indices)
    ignored = Tuple(new_axis(axis)::Int for axis in box.globally_ignored if !(axis in fixed_axes))
    entries = Dict(surviving)
    return CompoundBoundingBox{length(free)}(entries, selectors; selector_ignored=flags, ignored)
end

function _domain_pairs(values)
    values === nothing && return Pair[]
    values isa Pair && return [values]
    values isa NamedTuple || values isa AbstractDict ? collect(pairs(values)) : collect(values)
end

function _domain_fixed_values(fixed, count)
    result = Dict{Int,Real}()
    for pair in _domain_pairs(fixed)
        pair isa Pair || throw(ArgumentError("fixed domain values must be provided as pairs"))
        axis = pair.first
        axis isa Integer && !(axis isa Bool) || throw(ArgumentError("fixed domain axes must be integer indices"))
        1 <= axis <= count || throw(ArgumentError("fixed domain axis $axis must be between 1 and $count"))
        pair.second isa Real || throw(ArgumentError("fixed domain values must be real numbers"))
        haskey(result, axis) && throw(ArgumentError("fixed domain axis $axis was provided more than once"))
        result[Int(axis)] = pair.second
    end
    return result
end

function _sampling_interval(interval, step, center)
    step isa Real || throw(ArgumentError("grid step values must be real numbers"))
    step > 0 || throw(ArgumentError("grid step values must be positive"))
    lower, upper = interval
    if center
        lower = floor(lower + 0.5)
        upper = ceil(upper - 0.5)
    end
    count = ceil(Int, (upper - lower) / step) + 1
    return collect(range(lower; step, length=count))
end

function _grid_from_simple_box(box::BoundingBox{N}, fixed; step=1, center=true) where {N}
    fixed_values = _domain_fixed_values(fixed, N)
    for (axis, value) in fixed_values
        interval = bounds(box, axis)
        interval === nothing && continue
        interval[1] <= value <= interval[2] || throw(DomainError(value, "fixed input $axis lies outside the bounding box"))
    end
    sampled_axes = Tuple(axis for axis in bounded_axes(box) if !haskey(fixed_values, axis))
    missing = Tuple(axis for axis in ignored_axes(box) if !haskey(fixed_values, axis))
    isempty(missing) || throw(ArgumentError("fixed values are required for ignored axes $missing"))
    steps = step isa Real ? ntuple(_ -> step, length(sampled_axes)) : Tuple(step)
    length(steps) == length(sampled_axes) || throw(DimensionMismatch("grid step has length $(length(steps)), expected $(length(sampled_axes))"))
    axis_values = ntuple(index -> _sampling_interval(bounds(box, sampled_axes[index]), steps[index], center), length(sampled_axes))
    points = isempty(axis_values) ? ((),) : Iterators.product(axis_values...)
    T = promote_type(map(typeof, Base.values(fixed_values))..., map(eltype, axis_values)..., Float64)
    count = isempty(axis_values) ? 1 : prod(length, axis_values)
    result = Matrix{T}(undef, N, count)
    for (column, point) in enumerate(points)
        @inbounds for axis in 1:N
            fixed_index = get(fixed_values, axis, nothing)
            if fixed_index === nothing
                position = findfirst(==(axis), sampled_axes)::Int
                result[axis, column] = point[position]
            else
                result[axis, column] = fixed_index
            end
        end
    end
    return result
end

function grid_from_bounding_box(box::BoundingBox; step=1, center=true, selector=nothing, fixed=())
    selector === nothing || throw(ArgumentError("selector is only valid for a compound bounding box"))
    return _grid_from_simple_box(box, fixed; step, center)
end

function grid_from_bounding_box(box::CompoundBoundingBox{N,K}; step=1, center=true, selector=nothing, fixed=()) where {N,K}
    if selector === nothing
        length(box.boxes) == 1 || throw(ArgumentError("selector is required for a compound bounding box with multiple entries"))
        key, selected = only(box.boxes)
    else
        key = _compound_key(selector, K)
        selected = select_bounding_box(box, key)
    end
    fixed_values = _domain_fixed_values(fixed, N)
    for index in 1:K
        axis = box.selector_axes[index]
        value = key[index]
        haskey(fixed_values, axis) && !isequal(fixed_values[axis], value) && throw(ArgumentError("fixed selector axis $axis does not match selector value $value"))
        fixed_values[axis] = value
    end
    return _grid_from_simple_box(selected, fixed_values; step, center)
end
