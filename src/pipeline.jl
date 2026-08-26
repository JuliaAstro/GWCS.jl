"""Inclusive input-domain intervals for a WCS."""
struct BoundingBox{N,T<:Real}
    intervals::NTuple{N,Tuple{T,T}}
    order::Symbol
end

function BoundingBox(intervals::NTuple{N,Tuple{T,T}}; order=:F) where {N,T<:Real}
    all(interval -> interval[1] <= interval[2], intervals) || throw(ArgumentError("bounding-box lower limits must not exceed upper limits"))
    order in (:F, :C) || throw(ArgumentError("bounding-box order must be :F or :C"))
    canonical = order === :C ? reverse(intervals) : intervals
    return BoundingBox{N,T}(canonical, :F)
end

function BoundingBox(intervals::AbstractVector; order=:F)
    values = Tuple((promote(interval[1], interval[2])...,) for interval in intervals)
    T = promote_type(map(interval -> typeof(interval[1]), values)...)
    return BoundingBox(ntuple(i -> (T(values[i][1]), T(values[i][2])), length(values)); order)
end

@inline function Base.in(coords::StaticVector{N}, box::BoundingBox{N}) where {N}
    return all(i -> box.intervals[i][1] <= coords[i] <= box.intervals[i][2], 1:N)
end

"""A coordinate frame and the transform from that frame to the next step."""
struct Step{F<:AbstractCoordinateFrame,T}
    frame::F
    transform::T
end

Step(frame::AbstractCoordinateFrame) = Step(frame, nothing)

"""
    WCS(steps...; name="", pixel_shape=nothing, bounding_box=nothing)

Validated generalized WCS pipeline with precompiled forward and inverse roots.
"""
struct WCS{Nin,Nout,S,F,B,BB,PS}
    steps::S
    forward::F
    backward::B
    name::String
    pixel_shape::PS
    bounding_box::BB
end

function WCS(steps::Step...; name="", pixel_shape=nothing, bounding_box=nothing)
    length(steps) >= 2 || throw(ArgumentError("a WCS pipeline requires at least an input and output frame"))
    steps[end].transform === nothing || throw(ArgumentError("the final pipeline step must not have a transform"))
    all(step -> step.transform !== nothing, steps[1:(end - 1)]) || throw(ArgumentError("only the final pipeline step may omit its transform"))

    names = map(step -> step.frame.name, steps)
    length(unique(names)) == length(names) || throw(ArgumentError("frame names must be unique"))

    for i in 1:(length(steps) - 1)
        step = steps[i]
        next_step = steps[i + 1]
        input_dimension(step.transform) == frame_dimension(step.frame) || throw(DimensionMismatch("transform from frame $(step.frame.name) expects $(input_dimension(step.transform)) inputs but the frame has $(frame_dimension(step.frame)) axes"))
        output_dimension(step.transform) == frame_dimension(next_step.frame) || throw(DimensionMismatch("transform from frame $(step.frame.name) has $(output_dimension(step.transform)) outputs but frame $(next_step.frame.name) has $(frame_dimension(next_step.frame)) axes"))
    end

    compiled = ntuple(i -> i == length(steps) ? steps[i] : Step(steps[i].frame, _compile_transform_checked(steps[i].transform)), length(steps))
    forward = compose(map(step -> step.transform, compiled[1:(end - 1)])...)
    backward = _inverse_transform_checked(forward)
    Nin = frame_dimension(first(compiled).frame)
    Nout = frame_dimension(last(compiled).frame)
    bounding_box === nothing || length(bounding_box.intervals) == Nin || throw(DimensionMismatch("bounding box has $(length(bounding_box.intervals)) axes, expected $Nin"))
    pixel_shape === nothing || length(pixel_shape) == Nin || throw(DimensionMismatch("pixel shape has length $(length(pixel_shape)), expected $Nin"))
    shape = pixel_shape === nothing ? nothing : Tuple(Int.(pixel_shape))
    return WCS{Nin,Nout,typeof(compiled),typeof(forward),typeof(backward),typeof(bounding_box),typeof(shape)}(compiled, forward, backward, String(name), shape, bounding_box)
end

function fix_inputs(wcs::WCS{Nin}, fixed; atol=nothing, rtol=nothing) where {Nin}
    frame = input_frame(wcs)
    input_pairs = _fixed_pairs(fixed)
    indexed = map(input_pairs) do pair
        _fixed_index(pair.first, frame.axes_names, Nin) => pair.second
    end
    fixed_transform = fix_inputs(first(wcs.steps).transform, indexed; atol, rtol)

    for (index, value) in zip(filter(index -> fixed_transform.sources[index] < 0, 1:Nin), fixed_transform.fixed_values)
        if wcs.bounding_box !== nothing
            lower, upper = wcs.bounding_box.intervals[index]
            lower <= value <= upper || throw(DomainError(value, "fixed input $(frame.axes_names[index]) lies outside the WCS bounding box"))
        end
    end

    free_indices = fixed_transform.free_indices
    reduced_frame = _subset_frame(frame, free_indices)
    reduced_box = if wcs.bounding_box === nothing || isempty(free_indices)
        nothing
    else
        BoundingBox(map(index -> wcs.bounding_box.intervals[index], free_indices))
    end
    reduced_shape = wcs.pixel_shape === nothing ? nothing : map(index -> wcs.pixel_shape[index], free_indices)
    steps = (Step(reduced_frame, fixed_transform), Base.tail(wcs.steps)...)
    return WCS(steps...; name=wcs.name, pixel_shape=reduced_shape, bounding_box=reduced_box)
end

input_frame(wcs::WCS) = first(wcs.steps).frame
output_frame(wcs::WCS) = last(wcs.steps).frame
available_frames(wcs::WCS) = map(step -> step.frame.name, wcs.steps)

function _frame_index(wcs::WCS, frame::AbstractString)
    index = findfirst(==(frame), available_frames(wcs))
    index === nothing && throw(ArgumentError("frame $frame is not in the pipeline"))
    return index
end

_frame_index(wcs::WCS, frame::AbstractCoordinateFrame) = _frame_index(wcs, frame.name)

function get_transform(wcs::WCS, from_frame, to_frame)
    from_index = _frame_index(wcs, from_frame)
    to_index = _frame_index(wcs, to_frame)
    from_index == to_index && return IdentityTransform{frame_dimension(wcs.steps[from_index].frame)}()

    if from_index < to_index
        return compose(map(step -> step.transform, wcs.steps[from_index:(to_index - 1)])...)
    end

    transforms = map(step -> _inverse_transform_checked(step.transform), reverse(wcs.steps[to_index:(from_index - 1)]))
    any(transform -> transform isa NoInverse, transforms) && return NoInverse()
    return compose(transforms...)
end

function transform(wcs::WCS, from_frame, to_frame, coords::StaticVector)
    compiled = get_transform(wcs, from_frame, to_frame)
    compiled isa NoInverse && throw(ArgumentError("no compiled transform is available from $from_frame to $to_frame"))
    return compiled(coords)
end

function transform(wcs::WCS, from_frame, to_frame, coords::Real...)
    return transform(wcs, from_frame, to_frame, SVector(coords))
end

transform(wcs::WCS, from_frame, to_frame) = transform(wcs, from_frame, to_frame, SVector{0,Float64}())
