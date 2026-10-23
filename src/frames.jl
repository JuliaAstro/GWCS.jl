"""Abstract coordinate frame with `N` axes."""
abstract type AbstractCoordinateFrame{N} end

@enum AxisType::UInt8 begin
    SPATIAL
    SPECTRAL
    TIME
    STOKES
    PIXEL
    UNDEFINED
end

frame_dimension(::AbstractCoordinateFrame{N}) where {N} = N
axis_names(frame::AbstractCoordinateFrame) = frame.axes_names
axis_physical_types(frame::AbstractCoordinateFrame) = frame.axis_physical_types
axis_units(frame::AbstractCoordinateFrame) = frame.unit
axis_types(frame::AbstractCoordinateFrame) = frame.axis_types

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

"""Generic coordinate frame metadata."""
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

A frame with `N` axes known only by its name, like Python GWCS's `EmptyFrame`, which a step frame written as a bare
string becomes. It has no axis names, physical types, units, or types, and is written back as its name.
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

function Frame2D(; name, axes_names=("x", "y"), axes_order=(1, 2), axis_types=nothing, axis_physical_types=nothing, unit=("", ""), reference_frame=nothing)
    types = axis_types === nothing ? (SPATIAL, SPATIAL) : axis_types
    return CoordinateFrame{2}(; name, axes_names, axes_order, axis_types=types, axis_physical_types, unit, reference_frame)
end

"""Celestial coordinate frame metadata."""
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

"""One-dimensional spectral coordinate frame metadata."""
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

"""One-dimensional temporal coordinate frame metadata."""
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

"""One-dimensional Stokes coordinate frame metadata."""
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

"""Mixed coordinate frame composed from heterogeneous constituent frames."""
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

_identity_axes_order(order) = all(index -> order[index] == index, eachindex(order))

function _frame_to_transform(frame::AbstractCoordinateFrame{N}) where {N}
    order = frame.axes_order
    _identity_axes_order(order) && return IdentityTransform{N}()
    mapping = ntuple(position -> findfirst(==(position), order)::Int, N)
    return RemapAxes{N}(mapping)
end

function _transform_to_frame(frame::AbstractCoordinateFrame{N}) where {N}
    _identity_axes_order(frame.axes_order) && return IdentityTransform{N}()
    return RemapAxes{N}(frame.axes_order)
end

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
