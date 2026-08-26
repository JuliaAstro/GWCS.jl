"""Abstract coordinate frame with `N` axes."""
abstract type AbstractCoordinateFrame{N} end

frame_dimension(::AbstractCoordinateFrame{N}) where {N} = N

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
    return CoordinateFrame{N,typeof(reference_frame)}(String(name), Tuple(String.(axes_names)), Tuple(Int.(axes_order)), Tuple(String.(axis_physical_types)), Tuple(String.(unit)), reference_frame)
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
    return CelestialFrame{N,typeof(reference_frame)}(String(name), Tuple(String.(axes_names)), Tuple(Int.(axes_order)), Tuple(String.(axis_physical_types)), Tuple(String.(unit)), reference_frame)
end

CelestialFrame(; kwargs...) = CelestialFrame{2}(; kwargs...)
