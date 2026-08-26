"""Abstract coordinate frame with `N` axes."""
abstract type AbstractCoordinateFrame{N} end

frame_dimension(::AbstractCoordinateFrame{N}) where {N} = N

function _axes_order(values, count)
    order = Tuple(values)
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
    axis_physical_types::NTuple{N,String}
    unit::NTuple{N,String}
    reference_frame::R
end

function CoordinateFrame{N}(; name, axes_names=ntuple(i -> "x$(i - 1)", N), axes_order=ntuple(identity, N), axis_physical_types=ntuple(i -> "custom:$(axes_names[i])", N), unit=ntuple(i -> "", N), reference_frame=nothing) where {N}
    return CoordinateFrame{N,typeof(reference_frame)}(String(name), Tuple(String.(axes_names)), _axes_order(axes_order, N), Tuple(String.(axis_physical_types)), Tuple(String.(unit)), reference_frame)
end

Frame2D(; kwargs...) = CoordinateFrame{2}(; kwargs...)

"""Celestial coordinate frame metadata."""
struct CelestialFrame{N,R} <: AbstractCoordinateFrame{N}
    name::String
    axes_names::NTuple{N,String}
    axes_order::NTuple{N,Int}
    axis_physical_types::NTuple{N,String}
    unit::NTuple{N,String}
    reference_frame::R
end

function CelestialFrame{N}(; name, axes_names=ntuple(i -> i == 1 ? "lon" : "lat", N), axes_order=ntuple(identity, N), axis_physical_types=ntuple(i -> i == 1 ? "pos.eq.ra" : "pos.eq.dec", N), unit=ntuple(i -> "deg", N), reference_frame=:ICRS) where {N}
    return CelestialFrame{N,typeof(reference_frame)}(String(name), Tuple(String.(axes_names)), _axes_order(axes_order, N), Tuple(String.(axis_physical_types)), Tuple(String.(unit)), reference_frame)
end

CelestialFrame(; kwargs...) = CelestialFrame{2}(; kwargs...)

function _validate_frame_axes_order(frame::AbstractCoordinateFrame{N}) where {N}
    all(position -> position in frame.axes_order, 1:N) || throw(ArgumentError("frame $(frame.name) axes_order must be a permutation of 1:$N"))
    return frame
end

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

function _subset_frame(frame::CoordinateFrame, indices::NTuple{N,Int}) where {N}
    return CoordinateFrame{N}(; name=frame.name, axes_names=map(index -> frame.axes_names[index], indices),
        axes_order=_subset_axes_order(frame.axes_order, indices),
        axis_physical_types=map(index -> frame.axis_physical_types[index], indices),
        unit=map(index -> frame.unit[index], indices), reference_frame=frame.reference_frame)
end

function _subset_frame(frame::CelestialFrame, indices::NTuple{N,Int}) where {N}
    return CelestialFrame{N}(; name=frame.name, axes_names=map(index -> frame.axes_names[index], indices),
        axes_order=_subset_axes_order(frame.axes_order, indices),
        axis_physical_types=map(index -> frame.axis_physical_types[index], indices),
        unit=map(index -> frame.unit[index], indices), reference_frame=frame.reference_frame)
end
