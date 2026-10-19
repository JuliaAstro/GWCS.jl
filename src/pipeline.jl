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
struct WCS{Nin,Nout,S,F,B,BB,PS,IS}
    steps::S
    forward::F
    backward::B
    name::String
    pixel_shape::PS
    bounding_box::BB
    inverse_split::IS
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
    kernel = compose(map(step -> _numeric_transform(step.transform), compiled[1:(end - 1)])...)
    forward = compose(_frame_to_transform(first(compiled).frame), kernel, _transform_to_frame(last(compiled).frame))
    backward = _inverse_transform_checked(forward)
    Nin = frame_dimension(first(compiled).frame)
    Nout = frame_dimension(last(compiled).frame)
    inverse_split = backward isa NoInverse && Nin == Nout ? _numerical_inverse_split(forward) : nothing
    if bounding_box !== nothing
        bounding_box isa AbstractBoundingBox || throw(ArgumentError("bounding_box must be an AbstractBoundingBox or nothing"))
        input_dimension(bounding_box) == Nin || throw(DimensionMismatch("bounding box has input dimension $(input_dimension(bounding_box)), expected $Nin"))
    end
    pixel_shape === nothing || length(pixel_shape) == Nin || throw(DimensionMismatch("pixel shape has length $(length(pixel_shape)), expected $Nin"))
    shape = pixel_shape === nothing ? nothing : Tuple(map(pixel_shape) do value
        value === nothing ? nothing : Int(value)
    end)
    return WCS{Nin,Nout,typeof(compiled),typeof(forward),typeof(backward),typeof(bounding_box),typeof(shape),typeof(inverse_split)}(compiled, forward, backward, String(name), shape, bounding_box, inverse_split)
end

"""
    WCSBuilder(wcs)
    WCSBuilder(steps...; name="", pixel_shape=nothing, bounding_box=nothing)

Checked construction-time wrapper whose state is always a valid immutable
`WCS`. Use `WCS(builder)` to obtain the current compiled snapshot.
"""
mutable struct WCSBuilder
    _wcs::WCS
end

WCSBuilder(steps::Step...; name="", pixel_shape=nothing, bounding_box=nothing) = WCSBuilder(WCS(steps...; name, pixel_shape, bounding_box))
WCS(builder::WCSBuilder) = builder._wcs

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

_frame_name(frame::AbstractString) = String(frame)
_frame_name(frame::Symbol) = String(frame)
_frame_name(frame::AbstractCoordinateFrame) = frame.name
_frame_name(frame) = throw(ArgumentError("frame must be a name, symbol, or coordinate frame"))

function _find_frame_index(wcs::WCS, frame)
    name = _frame_name(frame)
    return findfirst(==(name), available_frames(wcs))
end

function _frame_index(wcs::WCS, frame)
    name = _frame_name(frame)
    index = _find_frame_index(wcs, name)
    index === nothing && throw(ArgumentError("frame $frame is not in the pipeline"))
    return index
end

"""Return the canonical frame stored in `wcs` for a name, symbol, or frame object."""
get_frame(wcs::WCS, frame) = wcs.steps[_frame_index(wcs, frame)].frame

_replace_step(steps::Tuple, index, step) = ntuple(i -> i == index ? step : steps[i], length(steps))
_insert_step(steps::Tuple, index, step) = (steps[1:(index - 1)]..., step, steps[index:end]...)
_remove_step(steps::Tuple, index) = (steps[1:(index - 1)]..., steps[(index + 1):end]...)

function _rebuild_wcs(wcs::WCS, steps::Tuple; pixel_shape=wcs.pixel_shape, bounding_box=wcs.bounding_box)
    return WCS(steps...; name=wcs.name, pixel_shape, bounding_box)
end

"""Return a new WCS with the transform between adjacent frames replaced."""
function set_transform(wcs::WCS, from_frame, to_frame, transform::AbstractTransform)
    from_index = _frame_index(wcs, from_frame)
    to_index = _frame_index(wcs, to_frame)
    from_index + 1 == to_index || throw(ArgumentError("frames $(_frame_name(from_frame)) and $(_frame_name(to_frame)) are not adjacent in forward order"))
    steps = _replace_step(wcs.steps, from_index, Step(wcs.steps[from_index].frame, transform))
    return _rebuild_wcs(wcs, steps)
end

"""
    insert_transform(wcs, frame, transform; after=false)

Return a new WCS with `transform` composed onto the edge before or after
`frame`. Frame metadata is not changed when coordinate values change.
"""
function insert_transform(wcs::WCS, frame, transform::AbstractTransform; after=false)
    index = _frame_index(wcs, frame)
    edge = after ? index : index - 1
    edge >= 1 || throw(ArgumentError("a transform cannot be inserted before the first frame"))
    edge < length(wcs.steps) || throw(ArgumentError("a transform cannot be inserted after the last frame"))
    current = wcs.steps[edge].transform
    combined = after ? compose(transform, current) : compose(current, transform)
    steps = _replace_step(wcs.steps, edge, Step(wcs.steps[edge].frame, combined))
    return _rebuild_wcs(wcs, steps)
end

prepend_transform(wcs::WCS, transform::AbstractTransform) = insert_transform(wcs, input_frame(wcs), transform; after=true)
append_transform(wcs::WCS, transform::AbstractTransform) = insert_transform(wcs, output_frame(wcs), transform)

struct _UnspecifiedPipelineMetadata end
const _UNSPECIFIED_PIPELINE_METADATA = _UnspecifiedPipelineMetadata()

function _prepended_metadata(existing, replacement, property)
    replacement !== _UNSPECIFIED_PIPELINE_METADATA && return replacement
    existing === nothing && return nothing
    throw(ArgumentError("prepending an input frame requires an explicit $property replacement; pass $property=nothing to clear it"))
end

"""
    insert_frame(wcs, input_frame, transform, output_frame)

Insert a new frame connected by `transform`. Exactly one endpoint must already
exist. The added transform changes the overall mapping unless it is an identity.
"""
function insert_frame(wcs::WCS, input, transform::AbstractTransform, output; pixel_shape=_UNSPECIFIED_PIPELINE_METADATA, bounding_box=_UNSPECIFIED_PIPELINE_METADATA)
    input_index = _find_frame_index(wcs, input)
    output_index = _find_frame_index(wcs, output)
    input_index === nothing || output_index === nothing || throw(ArgumentError("exactly one inserted-frame endpoint must already exist in the pipeline"))
    input_index === nothing && output_index === nothing && throw(ArgumentError("exactly one inserted-frame endpoint must already exist in the pipeline"))

    if input_index === nothing
        input isa AbstractCoordinateFrame || throw(ArgumentError("new input endpoint must be a coordinate frame"))
        steps = _insert_step(wcs.steps, output_index, Step(input, transform))
        if output_index == 1
            shape = _prepended_metadata(wcs.pixel_shape, pixel_shape, "pixel_shape")
            box = _prepended_metadata(wcs.bounding_box, bounding_box, "bounding_box")
            return _rebuild_wcs(wcs, steps; pixel_shape=shape, bounding_box=box)
        end
        pixel_shape === _UNSPECIFIED_PIPELINE_METADATA || throw(ArgumentError("pixel_shape can only be changed when prepending an input frame"))
        bounding_box === _UNSPECIFIED_PIPELINE_METADATA || throw(ArgumentError("bounding_box can only be changed when prepending an input frame"))
        return _rebuild_wcs(wcs, steps)
    end

    output isa AbstractCoordinateFrame || throw(ArgumentError("new output endpoint must be a coordinate frame"))
    pixel_shape === _UNSPECIFIED_PIPELINE_METADATA || throw(ArgumentError("pixel_shape can only be changed when prepending an input frame"))
    bounding_box === _UNSPECIFIED_PIPELINE_METADATA || throw(ArgumentError("bounding_box can only be changed when prepending an input frame"))
    current = wcs.steps[input_index].transform
    replaced = _replace_step(wcs.steps, input_index, Step(wcs.steps[input_index].frame, transform))
    step = current === nothing ? Step(output) : Step(output, current)
    steps = _insert_step(replaced, input_index + 1, step)
    return _rebuild_wcs(wcs, steps)
end

"""Remove an intermediate frame and compose its adjacent transforms."""
function remove_frame(wcs::WCS, frame)
    index = _frame_index(wcs, frame)
    1 < index < length(wcs.steps) || throw(ArgumentError("only an intermediate frame can be removed"))
    incoming = wcs.steps[index - 1].transform
    outgoing = wcs.steps[index].transform
    combined = compose(incoming, outgoing)
    replaced = _replace_step(wcs.steps, index - 1, Step(wcs.steps[index - 1].frame, combined))
    return _rebuild_wcs(wcs, _remove_step(replaced, index))
end

"""Return a new WCS with one frame renamed and all other frame metadata preserved."""
function rename_frame(wcs::WCS, frame, name)
    index = _frame_index(wcs, frame)
    new_name = String(name)
    isempty(new_name) && throw(ArgumentError("frame name must not be empty"))
    any(==(new_name), available_frames(wcs)) && new_name != wcs.steps[index].frame.name && throw(ArgumentError("frame name $new_name is already in the pipeline"))
    renamed = _rename_frame(wcs.steps[index].frame, new_name)
    steps = _replace_step(wcs.steps, index, Step(renamed, wcs.steps[index].transform))
    return _rebuild_wcs(wcs, steps)
end

function prepend_frame(wcs::WCS, frame::AbstractCoordinateFrame, transform::AbstractTransform; pixel_shape=_UNSPECIFIED_PIPELINE_METADATA, bounding_box=_UNSPECIFIED_PIPELINE_METADATA)
    return insert_frame(wcs, frame, transform, input_frame(wcs); pixel_shape, bounding_box)
end

append_frame(wcs::WCS, transform::AbstractTransform, frame::AbstractCoordinateFrame) = insert_frame(wcs, output_frame(wcs), transform, frame)

function _commit!(builder::WCSBuilder, operation, args...; kwargs...)
    candidate = operation(builder._wcs, args...; kwargs...)
    builder._wcs = candidate
    return builder
end

set_transform!(builder::WCSBuilder, args...) = _commit!(builder, set_transform, args...)
insert_transform!(builder::WCSBuilder, args...; kwargs...) = _commit!(builder, insert_transform, args...; kwargs...)
prepend_transform!(builder::WCSBuilder, args...) = _commit!(builder, prepend_transform, args...)
append_transform!(builder::WCSBuilder, args...) = _commit!(builder, append_transform, args...)
insert_frame!(builder::WCSBuilder, args...; kwargs...) = _commit!(builder, insert_frame, args...; kwargs...)
remove_frame!(builder::WCSBuilder, args...) = _commit!(builder, remove_frame, args...)
rename_frame!(builder::WCSBuilder, args...) = _commit!(builder, rename_frame, args...)
prepend_frame!(builder::WCSBuilder, args...; kwargs...) = _commit!(builder, prepend_frame, args...; kwargs...)
append_frame!(builder::WCSBuilder, args...) = _commit!(builder, append_frame, args...)

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
