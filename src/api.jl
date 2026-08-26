pixel_n_dim(::WCS{Nin}) where {Nin} = Nin
world_n_dim(::WCS{Nin,Nout}) where {Nin,Nout} = Nout

_coordinate_float_type(::Type{T}) where {T<:AbstractFloat} = T
_coordinate_float_type(::Type{T}) where {T<:Real} = Float64

function _coordinate_vector(coords::Tuple{Vararg{Real,N}}) where {N}
    T = _coordinate_float_type(promote_type(map(typeof, coords)...))
    return SVector{N,T}(coords)
end

_coordinate_vector(::Tuple{}) = SVector{0,Float64}()

function _filled_coordinate(::Val{N}, ::Type{T}, fill_value) where {N,T}
    return SVector{N,T}(ntuple(i -> T(fill_value), N))
end

@inline function _in_domain(wcs::WCS{N}, pixels::StaticVector{N}; with_bounding_box=true) where {N}
    all(isfinite, pixels) || return false
    return !with_bounding_box || wcs.bounding_box === nothing || pixels in wcs.bounding_box
end

"""
    in_domain(wcs, pixels...) -> Bool
    in_domain(wcs, pixels::AbstractMatrix) -> BitVector

Test whether finite pixel coordinates belong to the WCS input domain. Matrix columns are coordinates.
"""
function in_domain(wcs::WCS{N}, coords::AbstractVector) where {N}
    length(coords) == N || throw(DimensionMismatch("pixel has length $(length(coords)), expected $N"))
    T = _coordinate_float_type(eltype(coords))
    return _in_domain(wcs, SVector{N,T}(coords))
end

in_domain(wcs::WCS, coords::Tuple{Vararg{Real}}) = _in_domain(wcs, _coordinate_vector(coords))
in_domain(wcs::WCS, coords::Real...) = in_domain(wcs, coords)
in_domain(wcs::WCS{0}) = _in_domain(wcs, SVector{0,Float64}())

function in_domain(wcs::WCS{N}, pixels::AbstractMatrix) where {N}
    size(pixels, 1) == N || throw(DimensionMismatch("pixel matrix has $(size(pixels, 1)) rows, expected $N"))
    T = _coordinate_float_type(eltype(pixels))
    result = BitVector(undef, size(pixels, 2))
    @inbounds for (index, column) in enumerate(axes(pixels, 2))
        result[index] = _in_domain(wcs, SVector{N,T}(view(pixels, :, column)))
    end
    return result
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
pixel_to_world(wcs::WCS{0}; kwargs...) = _pixel_to_world(wcs, SVector{0,Float64}(); kwargs...)

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

@inline function _in_image(wcs::WCS{Nin,Nout}, world::StaticVector{Nout}; with_bounding_box=true, tolerance=1e-5, maxiter=30, initial_guess=nothing, max_backtracks=12) where {Nin,Nout}
    pixels = try
        _world_to_pixel(wcs, world; with_bounding_box=false, tolerance, maxiter, initial_guess, quiet=false, max_backtracks)
    catch error
        if error isa NoConvergence || error isa DomainError
            return false
        end
        rethrow()
    end
    return _in_domain(wcs, pixels; with_bounding_box)
end

"""
    in_image(wcs, world...; with_bounding_box=true, kwargs...) -> Bool
    in_image(wcs, worlds::AbstractMatrix; with_bounding_box=true, kwargs...) -> BitVector

Test whether world coordinates invert successfully to finite pixels in the WCS input domain. Matrix columns are coordinates.
"""
function in_image(wcs::WCS{Nin,Nout}, coords::AbstractVector; kwargs...) where {Nin,Nout}
    length(coords) == Nout || throw(DimensionMismatch("world coordinate has length $(length(coords)), expected $Nout"))
    T = _coordinate_float_type(eltype(coords))
    return _in_image(wcs, SVector{Nout,T}(coords); kwargs...)
end

in_image(wcs::WCS, coords::Tuple{Vararg{Real}}; kwargs...) = _in_image(wcs, _coordinate_vector(coords); kwargs...)
in_image(wcs::WCS, coords::Real...; kwargs...) = in_image(wcs, coords; kwargs...)
in_image(wcs::WCS{Nin,0}; kwargs...) where {Nin} = _in_image(wcs, SVector{0,Float64}(); kwargs...)

function in_image(wcs::WCS{Nin,Nout}, worlds::AbstractMatrix; kwargs...) where {Nin,Nout}
    size(worlds, 1) == Nout || throw(DimensionMismatch("world matrix has $(size(worlds, 1)) rows, expected $Nout"))
    T = _coordinate_float_type(eltype(worlds))
    result = BitVector(undef, size(worlds, 2))
    @inbounds for (index, column) in enumerate(axes(worlds, 2))
        result[index] = _in_image(wcs, SVector{Nout,T}(view(worlds, :, column)); kwargs...)
    end
    return result
end
