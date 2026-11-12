"""
    AbstractCoordinateFrame{N}

Supertype of every coordinate frame, with `N` axes. A frame names a step's coordinates and describes their axes:
[`axis_names`](@ref), [`axis_units`](@ref), [`axis_physical_types`](@ref), and [`axis_types`](@ref), each in transform
order.

Python equivalent: `gwcs.coordinate_frames.BaseCoordinateFrame`.
"""
abstract type AbstractCoordinateFrame{N} end

@enum AxisType::UInt8 begin
    SPATIAL
    SPECTRAL
    TIME
    STOKES
    PIXEL
    UNDEFINED
end

@doc """
    AxisType

The kind of a frame axis: `GWCS.SPATIAL`, `GWCS.SPECTRAL`, `GWCS.TIME`, `GWCS.STOKES`, `GWCS.PIXEL`, or
`GWCS.UNDEFINED`. Frames also accept these names as case-insensitive strings, with `"TEMPORAL"` for `TIME`, and keep any
other nonempty string as a custom type.

Python equivalent: the strings of `gwcs.coordinate_frames.CoordinateFrame.axes_type`.
""" AxisType

"[`AxisType`](@ref) of a spatial axis, such as a detector position or a sky coordinate."
SPATIAL
"[`AxisType`](@ref) of a spectral axis, such as wavelength or frequency."
SPECTRAL
"[`AxisType`](@ref) of a time axis."
TIME
"[`AxisType`](@ref) of a Stokes polarization axis."
STOKES
"[`AxisType`](@ref) of a pixel axis."
PIXEL
"[`AxisType`](@ref) of an axis of unspecified type, the default of a [`CoordinateFrame`](@ref)."
UNDEFINED

"""
    frame_dimension(frame) -> Int

Return the number of axes of `frame`.

Python equivalent: `gwcs.coordinate_frames.CoordinateFrame.naxes`.
"""
frame_dimension(::AbstractCoordinateFrame{N}) where {N} = N

# As in Python GWCS, a frame stores its axes in its own order (a CelestialFrame's is lon, lat), and `axes_order[i]` is
# the transform position of frame axis `i`. Coordinates and the axis accessors follow the transform order. The
# entries of a composite constituent's `axes_order` are positions in the whole composite, so ranks are used.
_axis_rank(order, index) = count(<(order[index]), order) + 1
_axis_at(order, rank) = findfirst(index -> _axis_rank(order, index) == rank, ntuple(identity, length(order)))::Int
_to_transform_order(values::NTuple{N,Any}, order::NTuple{N,Int}) where {N} = ntuple(rank -> values[_axis_at(order, rank)], Val(N))

"""
    axis_names(frame) -> NTuple{N,String}

Return the names of the frame's axes in transform order.

Python equivalent: `gwcs.coordinate_frames.CoordinateFrame.axes_names`.
"""
axis_names(frame::AbstractCoordinateFrame) = _to_transform_order(frame.axes_names, frame.axes_order)

"""
    axis_physical_types(frame) -> NTuple{N,String}

Return the physical types of the frame's axes in transform order, as IVOA UCD strings such as `"pos.eq.ra"` or
`"custom:<name>"`.

Python equivalent: `gwcs.coordinate_frames.CoordinateFrame.axis_physical_types`.
"""
axis_physical_types(frame::AbstractCoordinateFrame) = _to_transform_order(frame.axis_physical_types, frame.axes_order)

"""
    axis_units(frame) -> NTuple{N,String}

Return the units of the frame's axes in transform order, as strings such as `"deg"`, empty for a unitless axis.

Python equivalent: `gwcs.coordinate_frames.CoordinateFrame.unit`, which holds astropy units rather than strings.
"""
axis_units(frame::AbstractCoordinateFrame) = _to_transform_order(frame.unit, frame.axes_order)

"""
    axis_types(frame) -> Tuple

Return the types of the frame's axes in transform order. Each entry is an [`AxisType`](@ref), such as `GWCS.SPATIAL`,
or a string for a custom type.

Python equivalent: `gwcs.coordinate_frames.CoordinateFrame.axes_type`, which holds strings.
"""
axis_types(frame::AbstractCoordinateFrame) = _to_transform_order(frame.axis_types, frame.axes_order)

_frame_values(values) = values isa Tuple || values isa AbstractVector ? Tuple(values) : (values,)

function _axis_type(value)
    value isa AxisType && return value
    value isa AbstractString || throw(ArgumentError("axis type must be an AxisType or string"))
    isempty(value) && throw(ArgumentError("custom axis type must not be empty"))
    normalized = uppercase(value)
    normalized == "SPATIAL" && return SPATIAL
    normalized == "SPECTRAL" && return SPECTRAL
    normalized in ("TIME", "TEMPORAL") && return TIME
    normalized == "STOKES" && return STOKES
    normalized == "PIXEL" && return PIXEL
    normalized == "UNDEFINED" && return UNDEFINED
    return String(value)
end

function _axis_types(values, count)
    types = _frame_values(values)
    length(types) == count || throw(DimensionMismatch("axis_types has length $(length(types)), expected $count"))
    return map(_axis_type, types)
end

function _frame_strings(values, count, property)
    strings = Tuple(String.(_frame_values(values)))
    length(strings) == count || throw(DimensionMismatch("$property has length $(length(strings)), expected $count"))
    return strings
end

function _axes_order(values, count)
    order = _frame_values(values)
    length(order) == count || throw(DimensionMismatch("axes_order has length $(length(order)), expected $count"))
    all(value -> value isa Integer && !(value isa Bool), order) || throw(ArgumentError("axes_order entries must be integers"))
    result = Tuple(Int.(order))
    all(>(0), result) || throw(ArgumentError("axes_order entries must be positive"))
    length(unique(result)) == count || throw(ArgumentError("axes_order entries must be unique"))
    return result
end

"""
    CoordinateFrame{N}(; name, axes_names=("x0", "x1", ...), axes_order=(1, 2, ...), axis_types=(UNDEFINED, ...),
                       axis_physical_types=nothing, unit=("", ...), reference_frame=nothing)

A frame of `N` axes of any kind, such as pixels or focal-plane coordinates. `axes_names`, `axis_types`,
`axis_physical_types`, and `unit` each hold one entry per axis, in the frame's own order. Axis types are
[`AxisType`](@ref)s or custom strings, and physical types default to `"custom:<axis name>"`.

`axes_order[i]` is the transform position of frame axis `i`, so it labels coordinates without reordering them; the
accessors such as [`axis_names`](@ref) return transform order. `reference_frame` holds metadata that describes the
coordinates, kept for round trips.

Python equivalent: `gwcs.coordinate_frames.CoordinateFrame`.
"""
struct CoordinateFrame{N,R} <: AbstractCoordinateFrame{N}
    name::String
    axes_names::NTuple{N,String}
    axes_order::NTuple{N,Int}
    axis_types::NTuple{N,Union{AxisType,String}}
    axis_physical_types::NTuple{N,String}
    unit::NTuple{N,String}
    reference_frame::R
end

function CoordinateFrame{N,R}(name, axes_names, axes_order, axis_physical_types, unit, reference_frame) where {N,R}
    types = ntuple(_ -> UNDEFINED, N)
    return CoordinateFrame{N,R}(String(name), _frame_strings(axes_names, N, "axes_names"), _axes_order(axes_order, N), types,
        _frame_strings(axis_physical_types, N, "axis_physical_types"), _frame_strings(unit, N, "unit"), reference_frame)
end

CoordinateFrame(name, axes_names::NTuple{N}, axes_order, axis_physical_types, unit, reference_frame::R) where {N,R} = CoordinateFrame{N,R}(name, axes_names, axes_order, axis_physical_types, unit, reference_frame)

function CoordinateFrame{N}(; name, axes_names=ntuple(i -> "x$(i - 1)", N), axes_order=ntuple(identity, N), axis_types=ntuple(_ -> UNDEFINED, N), axis_physical_types=nothing, unit=ntuple(_ -> "", N), reference_frame=nothing) where {N}
    names = _frame_strings(axes_names, N, "axes_names")
    types = _axis_types(axis_types, N)
    physical_values = axis_physical_types === nothing ? ntuple(i -> isempty(names[i]) ? "custom:$(types[i])" : "custom:$(names[i])", N) : axis_physical_types
    physical_types = _frame_strings(physical_values, N, "axis_physical_types")
    units = _frame_strings(unit, N, "unit")
    return CoordinateFrame{N,typeof(reference_frame)}(String(name), names, _axes_order(axes_order, N), types, physical_types, units, reference_frame)
end

"""
    EmptyFrame{N}(name)

A frame with `N` axes known only by its name, which a step frame written as a bare string in metadata becomes. It has no
axis names, physical types, units, or types, and is written back as its name.

Python equivalent: `gwcs.coordinate_frames.EmptyFrame`.
"""
struct EmptyFrame{N} <: AbstractCoordinateFrame{N}
    name::String
    axes_names::NTuple{N,String}
    axes_order::NTuple{N,Int}
    axis_types::NTuple{N,AxisType}
    axis_physical_types::NTuple{N,String}
    unit::NTuple{N,String}
    EmptyFrame{N}(name) where {N} = new{N}(name, ntuple(_ -> "", N), ntuple(identity, N), ntuple(_ -> UNDEFINED, N), ntuple(_ -> "", N), ntuple(_ -> "", N))
end

"""
    Frame2D(; name, axes_names=("x", "y"), axes_order=(1, 2), axis_types=nothing, axis_physical_types=nothing,
            unit=("", ""), reference_frame=nothing) -> CoordinateFrame{2}

Return a two-axis [`CoordinateFrame`](@ref) whose axes are spatial unless `axis_types` says otherwise, such as a detector
or focal-plane frame. Metadata writes a two-axis spatial frame without a reference frame as `frame2d`.

Python equivalent: `gwcs.coordinate_frames.Frame2D`.
"""
function Frame2D(; name, axes_names=("x", "y"), axes_order=(1, 2), axis_types=nothing, axis_physical_types=nothing, unit=("", ""), reference_frame=nothing)
    types = axis_types === nothing ? (SPATIAL, SPATIAL) : axis_types
    return CoordinateFrame{2}(; name, axes_names, axes_order, axis_types=types, axis_physical_types, unit, reference_frame)
end

"""
    CelestialFrame(; name, axes_names=("lon", "lat"), axes_order=(1, 2), axis_physical_types=("pos.eq.ra", "pos.eq.dec"),
                   unit=("deg", "deg"), reference_frame=:ICRS)
    CelestialFrame{N}(; kwargs...)

A frame of sky positions: longitude and latitude in that frame order, and with `N = 3` a distance. Its axes are
spatial. `reference_frame` identifies the celestial system; frames read from FITS or ASDF metadata hold the system it
describes, such as FK5 with its equinox. Per-axis keywords are in frame order; [`CoordinateFrame`](@ref) describes
`axes_order`.

Python equivalent: `gwcs.coordinate_frames.CelestialFrame`.
"""
struct CelestialFrame{N,R} <: AbstractCoordinateFrame{N}
    name::String
    axes_names::NTuple{N,String}
    axes_order::NTuple{N,Int}
    axis_types::NTuple{N,Union{AxisType,String}}
    axis_physical_types::NTuple{N,String}
    unit::NTuple{N,String}
    reference_frame::R
end

function CelestialFrame{N,R}(name, axes_names, axes_order, axis_physical_types, unit, reference_frame) where {N,R}
    types = ntuple(_ -> SPATIAL, N)
    return CelestialFrame{N,R}(String(name), _frame_strings(axes_names, N, "axes_names"), _axes_order(axes_order, N), types,
        _frame_strings(axis_physical_types, N, "axis_physical_types"), _frame_strings(unit, N, "unit"), reference_frame)
end

CelestialFrame(name, axes_names::NTuple{N}, axes_order, axis_physical_types, unit, reference_frame::R) where {N,R} = CelestialFrame{N,R}(name, axes_names, axes_order, axis_physical_types, unit, reference_frame)

function CelestialFrame{N}(; name, axes_names=ntuple(i -> i == 1 ? "lon" : i == 2 ? "lat" : "distance", N), axes_order=ntuple(identity, N), axis_physical_types=nothing, unit=ntuple(i -> i <= 2 ? "deg" : "", N), reference_frame=:ICRS) where {N}
    names = _frame_strings(axes_names, N, "axes_names")
    types = ntuple(_ -> SPATIAL, N)
    physical_values = axis_physical_types === nothing ? ntuple(i -> i == 1 ? "pos.eq.ra" : i == 2 ? "pos.eq.dec" : "pos.distance", N) : axis_physical_types
    physical_types = _frame_strings(physical_values, N, "axis_physical_types")
    units = _frame_strings(unit, N, "unit")
    return CelestialFrame{N,typeof(reference_frame)}(String(name), names, _axes_order(axes_order, N), types, physical_types, units, reference_frame)
end

CelestialFrame(; kwargs...) = CelestialFrame{2}(; kwargs...)

function _spectral_physical_type(unit)
    normalized = lowercase(strip(String(unit)))
    normalized in ("hz", "khz", "mhz", "ghz", "thz") && return "em.freq"
    normalized in ("m", "cm", "mm", "um", "micron", "microns", "nm", "angstrom", "aa") && return "em.wl"
    normalized in ("j", "joule", "joules", "ev", "kev", "mev", "gev", "erg") && return "em.energy"
    normalized in ("m/s", "m s-1", "km/s", "km s-1") && return "spect.dopplerVeloc"
    return "custom:spectral"
end

"""
    SpectralFrame(; name="spectral", axes_names=("spectral",), axes_order=(1,), axis_physical_types=nothing, unit=("",),
                  reference_frame=nothing, reference_position="geocenter")

A one-axis frame of a spectral coordinate. Without `axis_physical_types`, the unit sets the physical type: `"em.wl"` for
a length, `"em.freq"` for a frequency, `"em.energy"` for an energy, `"spect.dopplerVeloc"` for a velocity, and
`"custom:spectral"` otherwise. `reference_position` is `"geocenter"`, `"barycenter"`, or `"heliocenter"`.
[`CoordinateFrame`](@ref) describes the shared keywords.

Python equivalent: `gwcs.coordinate_frames.SpectralFrame`.
"""
struct SpectralFrame{R} <: AbstractCoordinateFrame{1}
    name::String
    axes_names::NTuple{1,String}
    axes_order::NTuple{1,Int}
    axis_types::NTuple{1,AxisType}
    axis_physical_types::NTuple{1,String}
    unit::NTuple{1,String}
    reference_frame::R
    reference_position::String
end

function SpectralFrame(; name="spectral", axes_names=("spectral",), axes_order=(1,), axis_physical_types=nothing, unit=("",), reference_frame=nothing, reference_position="geocenter")
    names = _frame_strings(axes_names, 1, "axes_names")
    units = _frame_strings(unit isa AbstractString ? (unit,) : unit, 1, "unit")
    physical_types = axis_physical_types === nothing ? (_spectral_physical_type(only(units)),) : _frame_strings(axis_physical_types, 1, "axis_physical_types")
    position = String(reference_position)
    position in ("geocenter", "barycenter", "heliocenter") || throw(ArgumentError("unsupported spectral reference position $position"))
    return SpectralFrame{typeof(reference_frame)}(String(name), names, _axes_order(axes_order, 1), (SPECTRAL,), physical_types, units, reference_frame, position)
end

"""
    TemporalFrame(; name="temporal", axes_names=("time",), axes_order=(1,), axis_physical_types=("time",), unit=("s",),
                  reference_frame=nothing, reference_epoch=nothing, time_scale="", time_format="", location=nothing)

A one-axis frame whose coordinate is a time offset in `unit` from `reference_epoch`, which defaults to
`reference_frame`. [`offset_to_time`](@ref) and [`time_to_offset`](@ref) convert between offsets and times.
`time_scale` and `time_format` describe the epoch, for example `"utc"` and `"mjd"`, and `location` is the observer's
location for topocentric times. [`CoordinateFrame`](@ref) describes the shared keywords.

Python equivalent: `gwcs.coordinate_frames.TemporalFrame`.
"""
struct TemporalFrame{E,L,R} <: AbstractCoordinateFrame{1}
    name::String
    axes_names::NTuple{1,String}
    axes_order::NTuple{1,Int}
    axis_types::NTuple{1,AxisType}
    axis_physical_types::NTuple{1,String}
    unit::NTuple{1,String}
    reference_frame::R
    reference_epoch::E
    time_scale::String
    time_format::String
    location::L
end

function TemporalFrame(; name="temporal", axes_names=("time",), axes_order=(1,), axis_physical_types=("time",), unit=("s",), reference_frame=nothing, reference_epoch=nothing, time_scale="", time_format="", location=nothing)
    names = _frame_strings(axes_names, 1, "axes_names")
    units = _frame_strings(unit isa AbstractString ? (unit,) : unit, 1, "unit")
    physical_types = _frame_strings(axis_physical_types, 1, "axis_physical_types")
    epoch = reference_epoch === nothing ? reference_frame : reference_epoch
    return TemporalFrame{typeof(epoch),typeof(location),typeof(reference_frame)}(String(name), names, _axes_order(axes_order, 1), (TIME,), physical_types, units, reference_frame, epoch, String(time_scale), String(time_format), location)
end

function _temporal_days_per_unit(unit)
    normalized = lowercase(strip(String(unit)))
    normalized in ("d", "day", "days") && return 1.0
    normalized in ("h", "hr", "hour", "hours") && return 1 / 24
    normalized in ("min", "minute", "minutes") && return 1 / 1440
    normalized in ("s", "sec", "second", "seconds") && return 1 / 86400
    throw(ArgumentError("cannot convert temporal unit $unit to days without an optional time integration"))
end

"""
    offset_to_time(frame::TemporalFrame, offset)

Return the time `offset` after the frame's reference epoch. With a numeric epoch in the `"jd"` or `"mjd"` format, the
result is a Julian or modified Julian date, with `offset` converted to days from the frame's unit (`d`, `h`, `min`, or
`s`). Otherwise it is `reference_epoch + offset`. A frame without a reference epoch, or with another numeric format,
throws an `ArgumentError`.
"""
function offset_to_time(frame::TemporalFrame, offset)
    frame.reference_epoch === nothing && throw(ArgumentError("temporal frame $(frame.name) has no reference epoch"))
    if frame.reference_epoch isa Real && lowercase(frame.time_format) in ("jd", "mjd")
        return frame.reference_epoch + offset * _temporal_days_per_unit(only(frame.unit))
    elseif frame.reference_epoch isa Real && !isempty(frame.time_format)
        throw(ArgumentError("numeric temporal format $(frame.time_format) requires an optional time integration"))
    end
    applicable(+, frame.reference_epoch, offset) || throw(ArgumentError("temporal reference $(typeof(frame.reference_epoch)) does not support offset arithmetic"))
    return frame.reference_epoch + offset
end

"""
    time_to_offset(frame::TemporalFrame, time)

Return the offset of `time` from the frame's reference epoch, in the frame's unit; the inverse of
[`offset_to_time`](@ref).
"""
function time_to_offset(frame::TemporalFrame, value)
    frame.reference_epoch === nothing && throw(ArgumentError("temporal frame $(frame.name) has no reference epoch"))
    if frame.reference_epoch isa Real && lowercase(frame.time_format) in ("jd", "mjd")
        return (value - frame.reference_epoch) / _temporal_days_per_unit(only(frame.unit))
    elseif frame.reference_epoch isa Real && !isempty(frame.time_format)
        throw(ArgumentError("numeric temporal format $(frame.time_format) requires an optional time integration"))
    end
    applicable(-, value, frame.reference_epoch) || throw(ArgumentError("temporal reference $(typeof(frame.reference_epoch)) does not support offset arithmetic"))
    return value - frame.reference_epoch
end

"""
    StokesFrame(; name="stokes", axes_names=("stokes",), axes_order=(1,), axis_physical_types=("phys.polarization.stokes",))

A unitless one-axis frame whose coordinate selects a Stokes parameter, numbered as in FITS, such as 1 to 4 for I, Q,
U, and V. [`CoordinateFrame`](@ref) describes the shared keywords.

Python equivalent: `gwcs.coordinate_frames.StokesFrame`.
"""
struct StokesFrame <: AbstractCoordinateFrame{1}
    name::String
    axes_names::NTuple{1,String}
    axes_order::NTuple{1,Int}
    axis_types::NTuple{1,AxisType}
    axis_physical_types::NTuple{1,String}
    unit::NTuple{1,String}
    reference_frame::Nothing
end

function StokesFrame(; name="stokes", axes_names=("stokes",), axes_order=(1,), axis_physical_types=("phys.polarization.stokes",))
    names = _frame_strings(axes_names, 1, "axes_names")
    physical_types = _frame_strings(axis_physical_types, 1, "axis_physical_types")
    return StokesFrame(String(name), names, _axes_order(axes_order, 1), (STOKES,), physical_types, ("",), nothing)
end

"""
    CompositeFrame(frames; name="composite")

A frame combining other frames, such as a celestial and a spectral frame for an image cube. Its axes are the
constituents' axes, concatenated in the order given. Each constituent's `axes_order` gives its axes' transform
positions in the whole composite, and together these must be a permutation of `1:N`, where `N` is the total number of
axes. Nested composites flatten.

Python equivalent: `gwcs.coordinate_frames.CompositeFrame`.
"""
struct CompositeFrame{N,F} <: AbstractCoordinateFrame{N}
    name::String
    frames::F
    axes_names::NTuple{N,String}
    axes_order::NTuple{N,Int}
    axis_types::NTuple{N,Union{AxisType,String}}
    axis_physical_types::NTuple{N,String}
    unit::NTuple{N,String}
    reference_frame::Nothing
end

_composite_leaves(frame::CompositeFrame) = frame.frames
_composite_leaves(frame::AbstractCoordinateFrame) = (frame,)
_flatten_composite_frames(::Tuple{}) = ()
_flatten_composite_frames(frames::Tuple) = (_composite_leaves(first(frames))..., _flatten_composite_frames(Base.tail(frames))...)
_flatten_frame_property(accessor, ::Tuple{}) = ()
_flatten_frame_property(accessor, frames::Tuple) = (accessor(first(frames))..., _flatten_frame_property(accessor, Base.tail(frames))...)

function CompositeFrame(frames; name="composite")
    supplied = Tuple(frames)
    values = _flatten_composite_frames(supplied)
    isempty(values) && throw(ArgumentError("a composite frame requires at least one constituent frame"))
    all(frame -> frame isa AbstractCoordinateFrame, values) || throw(ArgumentError("composite frame constituents must be coordinate frames"))
    N = sum(frame_dimension, values)
    axes_names = _flatten_frame_property(axis_names, values)
    axes_order = _flatten_frame_property(frame -> frame.axes_order, values)
    types = _flatten_frame_property(axis_types, values)
    physical_types = _flatten_frame_property(axis_physical_types, values)
    units = _flatten_frame_property(axis_units, values)
    order = _axes_order(axes_order, N)
    all(position -> position in order, 1:N) || throw(ArgumentError("composite frame axes_order must be a permutation of 1:$N"))
    return CompositeFrame{N,typeof(values)}(String(name), values, axes_names, order, types, physical_types, units, nothing)
end

function _validate_frame_axes_order(frame::AbstractCoordinateFrame{N}) where {N}
    all(position -> position in frame.axes_order, 1:N) || throw(ArgumentError("frame $(frame.name) axes_order must be a permutation of 1:$N"))
    return frame
end

function _rename_frame(frame::CoordinateFrame{N}, name) where {N}
    return CoordinateFrame{N}(; name, axes_names=frame.axes_names, axes_order=frame.axes_order, axis_types=frame.axis_types,
        axis_physical_types=frame.axis_physical_types, unit=frame.unit, reference_frame=frame.reference_frame)
end

function _rename_frame(frame::CelestialFrame{N}, name) where {N}
    return CelestialFrame{N}(; name, axes_names=frame.axes_names, axes_order=frame.axes_order,
        axis_physical_types=frame.axis_physical_types, unit=frame.unit, reference_frame=frame.reference_frame)
end

function _rename_frame(frame::SpectralFrame, name)
    return SpectralFrame(; name, axes_names=frame.axes_names, axes_order=frame.axes_order, axis_physical_types=frame.axis_physical_types,
        unit=frame.unit, reference_frame=frame.reference_frame, reference_position=frame.reference_position)
end

function _rename_frame(frame::TemporalFrame, name)
    return TemporalFrame(; name, axes_names=frame.axes_names, axes_order=frame.axes_order, axis_physical_types=frame.axis_physical_types,
        unit=frame.unit, reference_frame=frame.reference_frame, reference_epoch=frame.reference_epoch, time_scale=frame.time_scale,
        time_format=frame.time_format, location=frame.location)
end

function _rename_frame(frame::StokesFrame, name)
    return StokesFrame(; name, axes_names=frame.axes_names, axes_order=frame.axes_order, axis_physical_types=frame.axis_physical_types)
end

_rename_frame(frame::CompositeFrame, name) = CompositeFrame(frame.frames; name)
_rename_frame(::EmptyFrame{N}, name) where {N} = EmptyFrame{N}(name)
_rename_frame(frame::AbstractCoordinateFrame, name) = throw(ArgumentError("renaming frame type $(typeof(frame)) is not supported"))

function _subset_axes_order(order, indices)
    selected = map(index -> order[index], indices)
    ranked = sort(collect(selected))
    return map(value -> findfirst(==(value), ranked)::Int, selected)
end

function _subset_frame(frame::CoordinateFrame, indices::NTuple{N,Int}, axes_order=_subset_axes_order(frame.axes_order, indices)) where {N}
    return CoordinateFrame{N}(; name=frame.name, axes_names=map(index -> frame.axes_names[index], indices), axes_order,
        axis_types=map(index -> frame.axis_types[index], indices), axis_physical_types=map(index -> frame.axis_physical_types[index], indices),
        unit=map(index -> frame.unit[index], indices), reference_frame=frame.reference_frame)
end

_subset_frame(frame::EmptyFrame, indices::NTuple{N,Int}, axes_order=indices) where {N} = EmptyFrame{N}(frame.name)

function _subset_frame(frame::CelestialFrame, indices::NTuple{N,Int}, axes_order=_subset_axes_order(frame.axes_order, indices)) where {N}
    return CelestialFrame{N}(; name=frame.name, axes_names=map(index -> frame.axes_names[index], indices), axes_order,
        axis_physical_types=map(index -> frame.axis_physical_types[index], indices), unit=map(index -> frame.unit[index], indices),
        reference_frame=frame.reference_frame)
end

_subset_frame(frame::SpectralFrame, ::Tuple{}, axes_order=()) = CoordinateFrame{0}(; name=frame.name)

function _subset_frame(frame::SpectralFrame, indices::NTuple{1,Int}, axes_order=_subset_axes_order(frame.axes_order, indices))
    return SpectralFrame(; name=frame.name, axes_names=frame.axes_names, axes_order, axis_physical_types=frame.axis_physical_types, unit=frame.unit, reference_frame=frame.reference_frame, reference_position=frame.reference_position)
end

_subset_frame(frame::TemporalFrame, ::Tuple{}, axes_order=()) = CoordinateFrame{0}(; name=frame.name)

function _subset_frame(frame::TemporalFrame, indices::NTuple{1,Int}, axes_order=_subset_axes_order(frame.axes_order, indices))
    return TemporalFrame(; name=frame.name, axes_names=frame.axes_names, axes_order, axis_physical_types=frame.axis_physical_types, unit=frame.unit,
        reference_frame=frame.reference_frame, reference_epoch=frame.reference_epoch, time_scale=frame.time_scale, time_format=frame.time_format, location=frame.location)
end

_subset_frame(frame::StokesFrame, ::Tuple{}, axes_order=()) = CoordinateFrame{0}(; name=frame.name)

function _subset_frame(frame::StokesFrame, indices::NTuple{1,Int}, axes_order=_subset_axes_order(frame.axes_order, indices))
    return StokesFrame(; name=frame.name, axes_names=frame.axes_names, axes_order, axis_physical_types=frame.axis_physical_types)
end

function _subset_composite_frame_dynamic(frame::CompositeFrame, indices)
    isempty(indices) && return CoordinateFrame{0}(; name=frame.name)
    selected_order = _subset_axes_order(frame.axes_order, indices)
    selected_positions = Dict(indices[index] => selected_order[index] for index in eachindex(indices))
    children = ()
    offset = 0
    for child in frame.frames
        child_indices = Tuple(index for index in 1:frame_dimension(child) if haskey(selected_positions, offset + index))
        if !isempty(child_indices)
            child_order = map(index -> selected_positions[offset + index], child_indices)
            children = (children..., _subset_frame(child, child_indices, child_order))
        end
        offset += frame_dimension(child)
    end
    return CompositeFrame(children; name=frame.name)
end

@generated function _subset_composite_frame(frame::CompositeFrame{N,F}, ::Val{indices}) where {N,F,indices}
    isempty(indices) && return :(CoordinateFrame{0}(; name=frame.name))
    child_types = F.parameters
    all(type -> type isa DataType && type <: AbstractCoordinateFrame, child_types) || return :(_subset_composite_frame_dynamic(frame, $(QuoteNode(indices))))

    children = Any[]
    offset = 0
    for (child_index, child_type) in enumerate(child_types)
        dimension = supertype(child_type).parameters[1]
        child_indices = Tuple(index for index in 1:dimension if offset + index in indices)
        if !isempty(child_indices)
            selected_positions = map(index -> findfirst(==(offset + index), indices)::Int, child_indices)
            child_order = Expr(:tuple, map(position -> :(selected_order[$position]), selected_positions)...)
            push!(children, :(_subset_frame(frame.frames[$child_index], $(QuoteNode(child_indices)), $child_order)))
        end
        offset += dimension
    end
    child_tuple = Expr(:tuple, children...)
    return quote
        selected_order = _subset_axes_order(frame.axes_order, $(QuoteNode(indices)))
        return CompositeFrame($child_tuple; name=frame.name)
    end
end

Base.@constprop :aggressive function _subset_frame(frame::CompositeFrame, indices::NTuple{N,Int}) where {N}
    issorted(indices) || throw(ArgumentError("composite frame subset indices must be sorted"))
    return _subset_composite_frame(frame, Val(indices))
end
