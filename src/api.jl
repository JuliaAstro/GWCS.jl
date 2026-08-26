pixel_n_dim(::WCS{Nin}) where {Nin} = Nin
world_n_dim(::WCS{Nin,Nout}) where {Nin,Nout} = Nout

_coordinate_float_type(::Type{T}) where {T<:AbstractFloat} = T
_coordinate_float_type(::Type{T}) where {T<:Real} = Float64

function _coordinate_vector(coords::Tuple{Vararg{Real,N}}) where {N}
    T = _coordinate_float_type(promote_type(map(typeof, coords)...))
    return SVector{N,T}(coords)
end

function _filled_coordinate(::Val{N}, ::Type{T}, fill_value) where {N,T}
    return SVector{N,T}(ntuple(i -> T(fill_value), N))
end

@inline function _pixel_to_world(wcs::WCS{Nin,Nout}, coords::StaticVector{Nin}; with_bounding_box=true, fill_value=NaN) where {Nin,Nout}
    T = _coordinate_float_type(eltype(coords))
    pixels = SVector{Nin,T}(coords)
    if with_bounding_box && wcs.bounding_box !== nothing && !(pixels in wcs.bounding_box)
        return _filled_coordinate(Val(Nout), T, fill_value)
    end
    return wcs.forward(pixels)
end

function pixel_to_world(wcs::WCS{Nin}, coords::AbstractVector; kwargs...) where {Nin}
    length(coords) == Nin || throw(DimensionMismatch("pixel has length $(length(coords)), expected $Nin"))
    T = _coordinate_float_type(eltype(coords))
    return _pixel_to_world(wcs, SVector{Nin,T}(coords); kwargs...)
end

pixel_to_world(wcs::WCS, coords::Tuple{Vararg{Real}}; kwargs...) = _pixel_to_world(wcs, _coordinate_vector(coords); kwargs...)
pixel_to_world(wcs::WCS, coords::Real...; kwargs...) = pixel_to_world(wcs, coords; kwargs...)

function pixel_to_world(wcs::WCS{Nin,Nout}, pixels::AbstractMatrix; kwargs...) where {Nin,Nout}
    size(pixels, 1) == Nin || throw(DimensionMismatch("pixel matrix has $(size(pixels, 1)) rows, expected $Nin"))
    size(pixels, 2) == 0 && return Matrix{_coordinate_float_type(eltype(pixels))}(undef, Nout, 0)
    first_result = pixel_to_world(wcs, view(pixels, :, 1); kwargs...)
    result = Matrix{eltype(first_result)}(undef, Nout, size(pixels, 2))
    result[:, 1] .= first_result
    @inbounds for j in 2:size(pixels, 2)
        result[:, j] .= pixel_to_world(wcs, view(pixels, :, j); kwargs...)
    end
    return result
end

@inline function _world_to_pixel(wcs::WCS{Nin,Nout}, coords::StaticVector{Nout}; with_bounding_box=true, fill_value=NaN, kwargs...) where {Nin,Nout}
    T = _coordinate_float_type(eltype(coords))
    world = SVector{Nout,T}(coords)
    wcs.backward isa NoInverse && return _numerical_inverse(wcs, world; with_bounding_box, fill_value, kwargs...)
    pixels = wcs.backward(world)
    if with_bounding_box && wcs.bounding_box !== nothing && !(pixels in wcs.bounding_box)
        return _filled_coordinate(Val(Nin), eltype(pixels), fill_value)
    end
    return pixels
end

function world_to_pixel(wcs::WCS{Nin,Nout}, coords::AbstractVector; kwargs...) where {Nin,Nout}
    length(coords) == Nout || throw(DimensionMismatch("world coordinate has length $(length(coords)), expected $Nout"))
    T = _coordinate_float_type(eltype(coords))
    return _world_to_pixel(wcs, SVector{Nout,T}(coords); kwargs...)
end

world_to_pixel(wcs::WCS, coords::Tuple{Vararg{Real}}; kwargs...) = _world_to_pixel(wcs, _coordinate_vector(coords); kwargs...)
world_to_pixel(wcs::WCS, coords::Real...; kwargs...) = world_to_pixel(wcs, coords; kwargs...)

function world_to_pixel(wcs::WCS{Nin,Nout}, worlds::AbstractMatrix; kwargs...) where {Nin,Nout}
    size(worlds, 1) == Nout || throw(DimensionMismatch("world matrix has $(size(worlds, 1)) rows, expected $Nout"))
    size(worlds, 2) == 0 && return Matrix{_coordinate_float_type(eltype(worlds))}(undef, Nin, 0)
    first_result = world_to_pixel(wcs, view(worlds, :, 1); kwargs...)
    result = Matrix{eltype(first_result)}(undef, Nin, size(worlds, 2))
    result[:, 1] .= first_result
    @inbounds for j in 2:size(worlds, 2)
        result[:, j] .= world_to_pixel(wcs, view(worlds, :, j); kwargs...)
    end
    return result
end
