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
    foreach(step -> _validate_frame_axes_order(step.frame), steps)

    for i in 1:(length(steps) - 1)
        step = steps[i]
        next_step = steps[i + 1]
        input_dimension(step.transform) == frame_dimension(step.frame) || throw(DimensionMismatch("transform from frame $(step.frame.name) expects $(input_dimension(step.transform)) inputs but the frame has $(frame_dimension(step.frame)) axes"))
        output_dimension(step.transform) == frame_dimension(next_step.frame) || throw(DimensionMismatch("transform from frame $(step.frame.name) has $(output_dimension(step.transform)) outputs but frame $(next_step.frame.name) has $(frame_dimension(next_step.frame)) axes"))
    end

    compiled = ntuple(i -> i == length(steps) ? steps[i] : Step(steps[i].frame, _compile_transform_checked(steps[i].transform)), length(steps))
    kernel = compose(map(step -> step.transform, compiled[1:(end - 1)])...)
    forward = compose(_frame_to_transform(first(compiled).frame), kernel, _transform_to_frame(last(compiled).frame))
    backward = _inverse_transform_checked(forward)
    Nin = frame_dimension(first(compiled).frame)
    Nout = frame_dimension(last(compiled).frame)
    if bounding_box !== nothing
        bounding_box isa AbstractBoundingBox || throw(ArgumentError("bounding_box must be an AbstractBoundingBox or nothing"))
        input_dimension(bounding_box) == Nin || throw(DimensionMismatch("bounding box has input dimension $(input_dimension(bounding_box)), expected $Nin"))
    end
    pixel_shape === nothing || length(pixel_shape) == Nin || throw(DimensionMismatch("pixel shape has length $(length(pixel_shape)), expected $Nin"))
    shape = pixel_shape === nothing ? nothing : Tuple(Int.(pixel_shape))
    return WCS{Nin,Nout,typeof(compiled),typeof(forward),typeof(backward),typeof(bounding_box),typeof(shape)}(compiled, forward, backward, String(name), shape, bounding_box)
end

function fix_inputs(wcs::WCS{Nin}, fixed; atol=nothing, rtol=nothing) where {Nin}
    frame = input_frame(wcs)
    input_pairs = _fixed_pairs(fixed)
    native_indexed = map(input_pairs) do pair
        _fixed_index(pair.first, frame.axes_names, Nin) => pair.second
    end
    transform_indexed = map(native_indexed) do pair
        frame.axes_order[pair.first] => pair.second
    end
    fixed_transform = fix_inputs(first(wcs.steps).transform, transform_indexed; atol, rtol)

    fixed_indices = first.(native_indexed)
    free_indices = Tuple(filter(index -> !(index in fixed_indices), 1:Nin))
    reduced_frame = _subset_frame(frame, free_indices)
    reduced_box = wcs.bounding_box === nothing ? nothing : _reduce_bounding_box(wcs.bounding_box, native_indexed)
    reduced_shape = wcs.pixel_shape === nothing ? nothing : map(index -> wcs.pixel_shape[index], free_indices)
    steps = (Step(reduced_frame, fixed_transform), Base.tail(wcs.steps)...)
    return WCS(steps...; name=wcs.name, pixel_shape=reduced_shape, bounding_box=reduced_box)
end

function with_bounding_box(wcs::WCS, bounding_box)
    return WCS(wcs.steps...; name=wcs.name, pixel_shape=wcs.pixel_shape, bounding_box)
end

function with_pixel_shape(wcs::WCS, pixel_shape)
    return WCS(wcs.steps...; name=wcs.name, pixel_shape, bounding_box=wcs.bounding_box)
end

function _named_domain_pairs(frame::AbstractCoordinateFrame{N}, fixed) where {N}
    return map(_domain_pairs(fixed)) do pair
        _fixed_index(pair.first, frame.axes_names, N) => pair.second
    end
end

function grid_from_bounding_box(wcs::WCS; step=1, center=true, selector=nothing, fixed=())
    wcs.bounding_box === nothing && throw(ArgumentError("WCS has no bounding box"))
    values = _named_domain_pairs(input_frame(wcs), fixed)
    return grid_from_bounding_box(wcs.bounding_box; step, center, selector, fixed=values)
end

input_frame(wcs::WCS) = first(wcs.steps).frame
output_frame(wcs::WCS) = last(wcs.steps).frame
available_frames(wcs::WCS) = map(step -> step.frame.name, wcs.steps)
axis_correlation_matrix(wcs::WCS) = axis_correlation_matrix(wcs.forward)
pixel_axis_names(wcs::WCS) = axis_names(input_frame(wcs))
world_axis_names(wcs::WCS) = axis_names(output_frame(wcs))
pixel_axis_physical_types(wcs::WCS) = axis_physical_types(input_frame(wcs))
world_axis_physical_types(wcs::WCS) = axis_physical_types(output_frame(wcs))
pixel_axis_units(wcs::WCS) = axis_units(input_frame(wcs))
world_axis_units(wcs::WCS) = axis_units(output_frame(wcs))
pixel_axis_types(wcs::WCS) = axis_types(input_frame(wcs))
world_axis_types(wcs::WCS) = axis_types(output_frame(wcs))

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

    kernel = if from_index < to_index
        compose(map(step -> step.transform, wcs.steps[from_index:(to_index - 1)])...)
    else
        transforms = map(step -> _inverse_transform_checked(step.transform), reverse(wcs.steps[to_index:(from_index - 1)]))
        any(transform -> transform isa NoInverse, transforms) && return NoInverse()
        compose(transforms...)
    end

    source = wcs.steps[from_index].frame
    destination = wcs.steps[to_index].frame
    return compose(_frame_to_transform(source), kernel, _transform_to_frame(destination))
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
