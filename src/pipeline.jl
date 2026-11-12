"""
    Step(frame, transform)
    Step(frame)

One step of a [`WCS`](@ref): a coordinate frame and the transform from its coordinates to the next step's frame. The
last step of a pipeline has no transform.

Python equivalent: `gwcs.wcs.Step`.
"""
struct Step{F<:AbstractCoordinateFrame,T}
    frame::F
    transform::T
    # Unmodeled metadata of the step and its frame, by path below the step: "", "frame", "frame/frames/2", and so on.
    pass_through::_PassThrough
end

Step(frame::AbstractCoordinateFrame, transform) = Step(frame, transform, _NO_PASS_THROUGH)
Step(frame::AbstractCoordinateFrame) = Step(frame, nothing)

# An edit that keeps a step's frame keeps the step's unmodeled metadata.
_with_transform(step::Step, transform) = Step(step.frame, transform, step.pass_through)

"""
    WCS(steps::Step...; name="", pixel_shape=nothing, bounding_box=nothing)

A WCS pipeline: a sequence of coordinate frames, each [`Step`](@ref) holding a frame and the transform to the next,
from the input (pixel) frame to the output (world) frame. Construction validates the pipeline and compiles its
transforms once, for fast repeated evaluation with [`pixel_to_world`](@ref) and [`world_to_pixel`](@ref). A `WCS` is
immutable; editing functions such as [`set_transform`](@ref) return a new one.

Each step's transform maps its frame's dimension to the next frame's. The last step has no transform, and frame names
are unique. An invalid pipeline throws an `ArgumentError` or a `DimensionMismatch`.

`pixel_shape` gives the image size along each pixel axis, with `nothing` for an axis without one. `bounding_box`, a
[`BoundingBox`](@ref) or [`CompoundBoundingBox`](@ref) over the pixel axes, limits the input domain.

When every step has an analytical inverse, the WCS inverse composes them. Otherwise `world_to_pixel` solves
numerically, which needs as many world as pixel axes.

# Examples
```jldoctest
julia> using GWCS

julia> using GWCS: Frame2D, Scale, Shift, Step

julia> wcs = WCS(Step(Frame2D(name="detector"), Shift(1.0) & Scale(2.0)), Step(Frame2D(name="world")));

julia> available_frames(wcs)
("detector", "world")

julia> pixel_n_dim(wcs), world_n_dim(wcs)
(2, 2)
```

Python equivalent: `gwcs.WCS`, which takes the steps as a list of `(frame, transform)` pairs.
"""
struct WCS{Nin,Nout,S,F,B,BB,PS,IS}
    steps::S
    forward::F
    backward::B
    name::String
    pixel_shape::PS
    bounding_box::BB
    inverse_split::IS
    pass_through::_PassThrough
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

    compiled = ntuple(i -> i == length(steps) ? steps[i] : _with_transform(steps[i], _compile_transform_checked(steps[i].transform)), length(steps))
    forward = compose(map(step -> _numeric_transform(step.transform), compiled[1:(end - 1)])...)
    backward = _inverse_transform_checked(forward)
    Nin = frame_dimension(first(compiled).frame)
    Nout = frame_dimension(last(compiled).frame)
    if bounding_box !== nothing
        bounding_box isa AbstractBoundingBox || throw(ArgumentError("bounding_box must be an AbstractBoundingBox or nothing"))
        input_dimension(bounding_box) == Nin || throw(DimensionMismatch("bounding box has input dimension $(input_dimension(bounding_box)), expected $Nin"))
    end
    pixel_shape === nothing || length(pixel_shape) == Nin || throw(DimensionMismatch("pixel shape has length $(length(pixel_shape)), expected $Nin"))
    shape = pixel_shape === nothing ? nothing : Tuple(map(pixel_shape) do value
        value === nothing ? nothing : Int(value)
    end)
    inverse_split = backward isa NoInverse && Nin == Nout ? _numerical_inverse_split(forward, bounding_box, shape) : nothing
    return WCS{Nin,Nout,typeof(compiled),typeof(forward),typeof(backward),typeof(bounding_box),typeof(shape),typeof(inverse_split)}(compiled, forward, backward, String(name), shape, bounding_box, inverse_split, _NO_PASS_THROUGH)
end

# Edits that keep the WCS keep its unmodeled metadata.
_with_pass_through(wcs::WCS, pass_through) = typeof(wcs)(wcs.steps, wcs.forward, wcs.backward, wcs.name, wcs.pixel_shape, wcs.bounding_box, wcs.inverse_split, pass_through)

function fix_inputs(wcs::WCS{Nin}, fixed; atol=nothing, rtol=nothing) where {Nin}
    frame = input_frame(wcs)
    indexed = map(pair -> _fixed_index(pair.first, axis_names(frame), Nin) => pair.second, _fixed_pairs(fixed))
    fixed_transform = fix_inputs(first(wcs.steps).transform, indexed; atol, rtol)

    fixed_indices = first.(indexed)
    free_indices = Tuple(filter(index -> !(index in fixed_indices), 1:Nin))
    fixed_frame_axes = map(index -> _axis_at(frame.axes_order, index), fixed_indices)
    reduced_frame = _subset_frame(frame, Tuple(filter(index -> !(index in fixed_frame_axes), 1:Nin)))
    reduced_box = wcs.bounding_box === nothing ? nothing : _reduce_bounding_box(wcs.bounding_box, indexed)
    reduced_shape = wcs.pixel_shape === nothing ? nothing : map(index -> wcs.pixel_shape[index], free_indices)
    steps = (Step(reduced_frame, fixed_transform), Base.tail(wcs.steps)...)
    return _with_pass_through(WCS(steps...; name=wcs.name, pixel_shape=reduced_shape, bounding_box=reduced_box), wcs.pass_through)
end

"""
    with_bounding_box(wcs, bounding_box) -> WCS

Return a copy of `wcs` whose input domain is `bounding_box`, a [`BoundingBox`](@ref) or [`CompoundBoundingBox`](@ref)
over the pixel axes, or that has no bounding box for `nothing`.

Python equivalent: setting `gwcs.WCS.bounding_box`, which changes the WCS in place.
"""
function with_bounding_box(wcs::WCS, bounding_box)
    return _with_pass_through(WCS(wcs.steps...; name=wcs.name, pixel_shape=wcs.pixel_shape, bounding_box), wcs.pass_through)
end

"""
    with_pixel_shape(wcs, pixel_shape) -> WCS

Return a copy of `wcs` with `pixel_shape`, the image size along each pixel axis in transform order, or with no pixel
shape for `nothing`. An entry of `nothing` leaves that axis without a size.

Python equivalent: setting `gwcs.WCS.pixel_shape`, which changes the WCS in place.
"""
function with_pixel_shape(wcs::WCS, pixel_shape)
    return _with_pass_through(WCS(wcs.steps...; name=wcs.name, pixel_shape, bounding_box=wcs.bounding_box), wcs.pass_through)
end

function _named_domain_pairs(frame::AbstractCoordinateFrame{N}, fixed) where {N}
    return map(_domain_pairs(fixed)) do pair
        _fixed_index(pair.first, axis_names(frame), N) => pair.second
    end
end

"""
    grid_from_bounding_box(wcs; step=1, center=true, selector=nothing, fixed=(), origin=0) -> Matrix
    grid_from_bounding_box(box; step=1, center=true, selector=nothing, fixed=()) -> Matrix

Return a regular grid of points covering a bounding box, one point per column, ready for [`pixel_to_world`](@ref). Each
bounded axis is sampled every `step`, a number or one per sampled axis. With `center=true`, the bounds are first rounded
to the nearest pixel centers, whole numbers; with `center=false`, sampling starts at the lower bound and the last point
can pass the upper one, as in Python GWCS.

`fixed` holds axis-value pairs for axes that are not sampled, and is required for an ignored axis. A compound bounding
box samples the entry for `selector`, which is required when it has more than one entry.

For a WCS, `fixed` can also name axes by the input frame's [`axis_names`](@ref), and `origin=1` makes `fixed` values and
the grid one-based; `selector` is a selector value, not a pixel. A WCS without a bounding box throws an `ArgumentError`.

Python equivalent: `gwcs.wcstools.grid_from_bounding_box`, which returns one array per axis rather than a matrix.
"""
function grid_from_bounding_box(wcs::WCS; step=1, center=true, selector=nothing, fixed=(), origin=0)
    wcs.bounding_box === nothing && throw(ArgumentError("WCS has no bounding box"))
    offset = _origin_offset(origin)
    values = map(pair -> pair.first => pair.second - offset, _named_domain_pairs(input_frame(wcs), fixed))
    grid = grid_from_bounding_box(wcs.bounding_box; step, center, selector, fixed=values)
    grid .+= offset
    return grid
end

"""
    input_frame(wcs) -> AbstractCoordinateFrame

Return the frame of the first step, whose coordinates are the WCS pixels.

Python equivalent: `gwcs.WCS.input_frame`.
"""
input_frame(wcs::WCS) = first(wcs.steps).frame
"""
    output_frame(wcs) -> AbstractCoordinateFrame

Return the frame of the last step, whose coordinates are the WCS world coordinates.

Python equivalent: `gwcs.WCS.output_frame`.
"""
output_frame(wcs::WCS) = last(wcs.steps).frame
"""
    available_frames(wcs) -> Tuple{Vararg{String}}

Return the names of the WCS frames, from the input frame to the output frame.

Python equivalent: `gwcs.WCS.available_frames`.
"""
available_frames(wcs::WCS) = map(step -> step.frame.name, wcs.steps)
"""
    axis_correlation_matrix(wcs) -> SMatrix{world_n_dim,pixel_n_dim,Bool}
    axis_correlation_matrix(transform) -> SMatrix{output_dimension,input_dimension,Bool}

Return which world axes depend on which pixel axes: entry `(i, j)` is `true` when world axis `i` can change with
pixel axis `j`, both in transform order. A `false` entry guarantees independence, and a transform that cannot tell
reports `true`. For a transform it is [`transform_dependencies`](@ref).

Python equivalent: `gwcs.WCS.axis_correlation_matrix`, and `astropy.modeling.separable.separability_matrix` for a
transform.
"""
axis_correlation_matrix(wcs::WCS) = axis_correlation_matrix(wcs.forward)
"""
    pixel_axis_names(wcs) -> NTuple{pixel_n_dim,String}

Return the input frame's [`axis_names`](@ref) in transform order.

Python equivalent: `gwcs.WCS.pixel_axis_names`.
"""
pixel_axis_names(wcs::WCS) = axis_names(input_frame(wcs))
"""
    world_axis_names(wcs) -> NTuple{world_n_dim,String}

Return the output frame's [`axis_names`](@ref) in transform order.

Python equivalent: `gwcs.WCS.world_axis_names`.
"""
world_axis_names(wcs::WCS) = axis_names(output_frame(wcs))
"""
    pixel_axis_physical_types(wcs) -> NTuple{pixel_n_dim,String}

Return the input frame's [`axis_physical_types`](@ref) in transform order.

Python equivalent: `gwcs.WCS.input_frame.axis_physical_types`.
"""
pixel_axis_physical_types(wcs::WCS) = axis_physical_types(input_frame(wcs))
"""
    world_axis_physical_types(wcs) -> NTuple{world_n_dim,String}

Return the output frame's [`axis_physical_types`](@ref) in transform order, as strings such as `"pos.eq.ra"`.

Python equivalent: `gwcs.WCS.world_axis_physical_types`.
"""
world_axis_physical_types(wcs::WCS) = axis_physical_types(output_frame(wcs))
"""
    pixel_axis_units(wcs) -> NTuple{pixel_n_dim,String}

Return the input frame's [`axis_units`](@ref) in transform order.

Python equivalent: `gwcs.WCS.input_frame.unit`, which holds astropy units rather than strings.
"""
pixel_axis_units(wcs::WCS) = axis_units(input_frame(wcs))
"""
    world_axis_units(wcs) -> NTuple{world_n_dim,String}

Return the output frame's [`axis_units`](@ref) in transform order, as strings such as `"deg"`.

Python equivalent: `gwcs.WCS.world_axis_units`.
"""
world_axis_units(wcs::WCS) = axis_units(output_frame(wcs))
"""
    pixel_axis_types(wcs) -> Tuple

Return the input frame's [`axis_types`](@ref) in transform order. Each entry is an [`AxisType`](@ref), such as
`GWCS.SPATIAL`, or a string for a custom type.

Python equivalent: `gwcs.WCS.input_frame.axes_type`, which holds strings.
"""
pixel_axis_types(wcs::WCS) = axis_types(input_frame(wcs))
"""
    world_axis_types(wcs) -> Tuple

Return the output frame's [`axis_types`](@ref) in transform order. Each entry is an [`AxisType`](@ref), such as
`GWCS.SPATIAL`, or a string for a custom type.

Python equivalent: `gwcs.WCS.output_frame.axes_type`, which holds strings.
"""
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

"""
    get_frame(wcs, frame) -> AbstractCoordinateFrame

Return the frame stored in `wcs` for `frame`, given as a name, a symbol, or a frame with that name. A frame not in the
pipeline throws an `ArgumentError`.
"""
get_frame(wcs::WCS, frame) = wcs.steps[_frame_index(wcs, frame)].frame

_replace_step(steps::Tuple, index, step) = ntuple(i -> i == index ? step : steps[i], length(steps))
_insert_step(steps::Tuple, index, step) = (steps[1:(index - 1)]..., step, steps[index:end]...)
_remove_step(steps::Tuple, index) = (steps[1:(index - 1)]..., steps[(index + 1):end]...)

function _rebuild_wcs(wcs::WCS, steps::Tuple; pixel_shape=wcs.pixel_shape, bounding_box=wcs.bounding_box)
    return _with_pass_through(WCS(steps...; name=wcs.name, pixel_shape, bounding_box), wcs.pass_through)
end

"""
    set_transform(wcs, from_frame, to_frame, transform) -> WCS

Return a copy of `wcs` whose transform from `from_frame` to the next frame, `to_frame`, is `transform`. Frames are given
by name or as frames, and frames that are not adjacent in that order throw an `ArgumentError`.

Python equivalent: `gwcs.WCS.set_transform`, which changes the WCS in place.
"""
function set_transform(wcs::WCS, from_frame, to_frame, transform::AbstractTransform)
    from_index = _frame_index(wcs, from_frame)
    to_index = _frame_index(wcs, to_frame)
    from_index + 1 == to_index || throw(ArgumentError("frames $(_frame_name(from_frame)) and $(_frame_name(to_frame)) are not adjacent in forward order"))
    steps = _replace_step(wcs.steps, from_index, _with_transform(wcs.steps[from_index], transform))
    return _rebuild_wcs(wcs, steps)
end

"""
    insert_transform(wcs, frame, transform; after=false) -> WCS

Return a copy of `wcs` with `transform` inserted at `frame`. With `after=false`, it runs last in the transform leading
into `frame`; with `after=true`, it runs first in the transform leaving `frame`. Frame metadata is unchanged, although
the coordinate values in the frames it affects change.

Python equivalent: `gwcs.WCS.insert_transform`, which changes the WCS in place.
"""
function insert_transform(wcs::WCS, frame, transform::AbstractTransform; after=false)
    index = _frame_index(wcs, frame)
    edge = after ? index : index - 1
    edge >= 1 || throw(ArgumentError("a transform cannot be inserted before the first frame"))
    edge < length(wcs.steps) || throw(ArgumentError("a transform cannot be inserted after the last frame"))
    current = wcs.steps[edge].transform
    combined = after ? compose(transform, current) : compose(current, transform)
    steps = _replace_step(wcs.steps, edge, _with_transform(wcs.steps[edge], combined))
    return _rebuild_wcs(wcs, steps)
end

"""
    prepend_transform(wcs, transform) -> WCS

Return a copy of `wcs` with `transform` applied to the pixels first, as `insert_transform(wcs, input_frame(wcs),
transform; after=true)` does.
"""
prepend_transform(wcs::WCS, transform::AbstractTransform) = insert_transform(wcs, input_frame(wcs), transform; after=true)

"""
    append_transform(wcs, transform) -> WCS

Return a copy of `wcs` with `transform` applied last, before the output frame, as
`insert_transform(wcs, output_frame(wcs), transform)` does.
"""
append_transform(wcs::WCS, transform::AbstractTransform) = insert_transform(wcs, output_frame(wcs), transform)

struct _UnspecifiedPipelineMetadata end
const _UNSPECIFIED_PIPELINE_METADATA = _UnspecifiedPipelineMetadata()

function _prepended_metadata(existing, replacement, property)
    replacement !== _UNSPECIFIED_PIPELINE_METADATA && return replacement
    existing === nothing && return nothing
    throw(ArgumentError("prepending an input frame requires an explicit $property replacement; pass $property=nothing to clear it"))
end

"""
    insert_frame(wcs, input, transform, output; pixel_shape, bounding_box) -> WCS

Return a copy of `wcs` with a new frame joined to the pipeline by `transform`. Exactly one of `input` and `output` is a
frame already in `wcs`, given by name or as a frame, and the other is the new frame.

When `output` is new, it is inserted right after `input`. `transform` connects them, and the transform that left
`input` now starts at `output`.

When `input` is new, it is inserted right before `output`. `transform` connects them, and the transform that reached
`output` now ends at `input`.

Either way the overall mapping changes unless `transform` is an identity. `pixel_shape` and `bounding_box` are accepted
only when the new frame becomes the WCS input frame. Then each value `wcs` has needs a replacement for the new
coordinates, or `nothing` to clear it.

Python equivalent: `gwcs.WCS.insert_frame`, which changes the WCS in place.
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
    replaced = _replace_step(wcs.steps, input_index, _with_transform(wcs.steps[input_index], transform))
    step = current === nothing ? Step(output) : Step(output, current)
    steps = _insert_step(replaced, input_index + 1, step)
    return _rebuild_wcs(wcs, steps)
end

"""
    remove_frame(wcs, frame) -> WCS

Return a copy of `wcs` without the intermediate `frame`, with the transforms into and out of it composed into one. The
overall mapping is unchanged. Removing the input or output frame throws an `ArgumentError`.
"""
function remove_frame(wcs::WCS, frame)
    index = _frame_index(wcs, frame)
    1 < index < length(wcs.steps) || throw(ArgumentError("only an intermediate frame can be removed"))
    incoming = wcs.steps[index - 1].transform
    outgoing = wcs.steps[index].transform
    combined = compose(incoming, outgoing)
    replaced = _replace_step(wcs.steps, index - 1, _with_transform(wcs.steps[index - 1], combined))
    return _rebuild_wcs(wcs, _remove_step(replaced, index))
end

"""
    rename_frame(wcs, frame, name) -> WCS

Return a copy of `wcs` with `frame` named `name` and otherwise unchanged. An empty name, or one another frame already
has, throws an `ArgumentError`.
"""
function rename_frame(wcs::WCS, frame, name)
    index = _frame_index(wcs, frame)
    new_name = String(name)
    isempty(new_name) && throw(ArgumentError("frame name must not be empty"))
    any(==(new_name), available_frames(wcs)) && new_name != wcs.steps[index].frame.name && throw(ArgumentError("frame name $new_name is already in the pipeline"))
    renamed = _rename_frame(wcs.steps[index].frame, new_name)
    steps = _replace_step(wcs.steps, index, Step(renamed, wcs.steps[index].transform, wcs.steps[index].pass_through))
    return _rebuild_wcs(wcs, steps)
end

"""
    prepend_frame(wcs, frame, transform; pixel_shape, bounding_box) -> WCS

Return a copy of `wcs` with `frame` as its new input frame and `transform` mapping it to the old one, as
`insert_frame(wcs, frame, transform, input_frame(wcs))` does. When `wcs` has a pixel shape or bounding box, pass their
replacements for the new pixel coordinates, or `nothing` to clear them.
"""
function prepend_frame(wcs::WCS, frame::AbstractCoordinateFrame, transform::AbstractTransform; pixel_shape=_UNSPECIFIED_PIPELINE_METADATA, bounding_box=_UNSPECIFIED_PIPELINE_METADATA)
    return insert_frame(wcs, frame, transform, input_frame(wcs); pixel_shape, bounding_box)
end

"""
    append_frame(wcs, transform, frame) -> WCS

Return a copy of `wcs` with `frame` as its new output frame and `transform` mapping the old output frame to it, as
`insert_frame(wcs, output_frame(wcs), transform, frame)` does.
"""
append_frame(wcs::WCS, transform::AbstractTransform, frame::AbstractCoordinateFrame) = insert_frame(wcs, output_frame(wcs), transform, frame)

"""
    get_transform(wcs, from_frame, to_frame) -> Union{AbstractTransform,NoInverse}

Return the transform between two frames of `wcs`, each given by name or as a frame. Toward the output frame it composes
the steps between them; back toward the input frame it composes their inverses, and gives [`NoInverse`](@ref) when a
step has none. The same frame twice gives an identity, and a frame not in the pipeline throws an `ArgumentError`.

Python equivalent: `gwcs.WCS.get_transform`.
"""
function get_transform(wcs::WCS, from_frame, to_frame)
    from_index = _frame_index(wcs, from_frame)
    to_index = _frame_index(wcs, to_frame)
    from_index == to_index && return IdentityTransform{frame_dimension(wcs.steps[from_index].frame)}()

    from_index < to_index && return compose(map(step -> step.transform, wcs.steps[from_index:(to_index - 1)])...)
    transforms = map(step -> _inverse_transform_checked(step.transform), reverse(wcs.steps[to_index:(from_index - 1)]))
    any(transform -> transform isa NoInverse, transforms) && return NoInverse()
    return compose(transforms...)
end

"""
    GWCS.transform(wcs, from_frame, to_frame, coords...) -> SVector

Evaluate [`get_transform`](@ref)`(wcs, from_frame, to_frame)` at a coordinate, given as separate numbers or a static
vector. Unlike [`pixel_to_world`](@ref), it applies no bounding box and no `fill_value`. A direction without an inverse
throws an `ArgumentError`.

Python equivalent: `gwcs.WCS.transform`, which applies the bounding box by default.
"""
function transform(wcs::WCS, from_frame, to_frame, coords::StaticVector)
    compiled = get_transform(wcs, from_frame, to_frame)
    compiled isa NoInverse && throw(ArgumentError("no compiled transform is available from $from_frame to $to_frame"))
    return compiled(coords)
end

function transform(wcs::WCS, from_frame, to_frame, coords::Real...)
    return transform(wcs, from_frame, to_frame, SVector(coords))
end

transform(wcs::WCS, from_frame, to_frame) = transform(wcs, from_frame, to_frame, SVector{0,Float64}())
