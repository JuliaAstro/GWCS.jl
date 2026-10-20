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

# `strict=false` returns `fill_value` for coordinates without a result (outside the bounding box or
# undefined under the transforms); `strict=true` throws for them instead. NaN input propagates as NaN.
@inline function _pixel_to_world(wcs::WCS{Nin,Nout}, coords::StaticVector{Nin}; with_bounding_box=true, fill_value=NaN, strict=false) where {Nin,Nout}
    T = _coordinate_float_type(eltype(coords))
    pixels = SVector{Nin,T}(coords)
    finite = all(isfinite, pixels)
    if with_bounding_box && wcs.bounding_box !== nothing && finite && !(pixels in wcs.bounding_box)
        strict && _throw_outside_box(pixels)
        return _filled_coordinate(Val(Nout), T, fill_value)
    end
    world = wcs.forward(pixels)
    # Transforms return NaN where a coordinate has no image; finite input that yields NaN is undefined.
    any(isnan, world) || return world
    strict || return _filled_coordinate(Val(Nout), T, fill_value)
    finite || return world
    _throw_undefined(wcs, pixels)
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

@inline function _world_to_pixel(wcs::WCS{Nin,Nout}, coords::StaticVector{Nout}; with_bounding_box=true, fill_value=NaN, strict=false, kwargs...) where {Nin,Nout}
    T = _coordinate_float_type(eltype(coords))
    world = SVector{Nout,T}(coords)
    pixels = wcs.backward isa NoInverse ? _numerical_solve(wcs, world; kwargs...) : wcs.backward(world)
    return _finish_inverse(wcs, world, pixels; with_bounding_box, fill_value, strict)
end

# Shared by the analytical and numerical inverses: undefined results are thrown or filled before the bounding box applies.
@inline function _finish_inverse(wcs::WCS{Nin}, world, pixels; with_bounding_box, fill_value, strict) where {Nin}
    if any(isnan, pixels)
        strict || return _filled_coordinate(Val(Nin), eltype(pixels), fill_value)
        all(isfinite, world) && _throw_undefined_inverse(wcs, world)
        return pixels
    end
    if with_bounding_box && wcs.bounding_box !== nothing && !(pixels in wcs.bounding_box)
        strict && _throw_outside_box_inverse(world, pixels)
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
        _world_to_pixel(wcs, world; with_bounding_box=false, fill_value=NaN, strict=false, tolerance, maxiter, initial_guess, quiet=false, max_backtracks)
    catch error
        # DomainError remains only for zero-dimensional inverses, which cannot carry NaN.
        error isa NoConvergence || error isa DomainError || rethrow()
        return false
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

@noinline _throw_outside_box(pixels) = throw(DomainError(Tuple(pixels), "pixel $(Tuple(pixels)) is outside the WCS bounding box"))
@noinline _throw_outside_box_inverse(world, pixels) = throw(DomainError(Tuple(world), "world coordinate $(Tuple(world)) maps to pixel $(Tuple(pixels)), outside the WCS bounding box"))

# Error path only: re-evaluate one coordinate step by step to report which transform is undefined there.
_describe_transform(transform) = string(nameof(typeof(transform)))
_describe_transform(transform::Union{IntermediateToNative,NativeToIntermediate}) = "$(nameof(typeof(transform)))($(nameof(typeof(transform.projection))))"

_undefined_leaf(transform, coords) = (transform, coords)
_undefined_leaf(transform::AnnotatedTransform, coords) = _undefined_leaf(transform.transform, coords)
_undefined_leaf(transform::WithInverse, coords) = _undefined_leaf(transform.forward, coords)

function _undefined_leaf(transform::ComposedTransform, coords)
    for child in transform.transforms
        next = child(coords)
        any(isnan, next) && return _undefined_leaf(child, coords)
        coords = next
    end
    return (transform, coords)
end

function _undefined_leaf(transform::ParallelTransform, coords)
    offset = 0
    for child in transform.transforms
        count = input_dimension(child)
        part = SVector{count}(ntuple(i -> coords[offset + i], count))
        any(isnan, child(part)) && return _undefined_leaf(child, part)
        offset += count
    end
    return (transform, coords)
end

function _undefined_message(subject, step_description, transform, coords)
    leaf, leaf_coords = _undefined_leaf(transform, coords)
    return "$subject: $(_describe_transform(leaf)) in $step_description is undefined at $(Tuple(leaf_coords))"
end

@noinline function _throw_undefined(wcs::WCS, pixels)
    subject = "pixel $(Tuple(pixels)) has no world coordinate"
    coords = _frame_to_transform(first(wcs.steps).frame)(pixels)
    for index in 1:(length(wcs.steps) - 1)
        step = wcs.steps[index]
        next = step.transform(coords)
        if any(isnan, next)
            description = "the step from frame \"$(step.frame.name)\" to \"$(wcs.steps[index + 1].frame.name)\""
            throw(DomainError(Tuple(pixels), _undefined_message(subject, description, step.transform, coords)))
        end
        coords = next
    end
    throw(DomainError(Tuple(pixels), subject))
end

@noinline function _throw_undefined_inverse(wcs::WCS, world)
    subject = "world coordinate $(Tuple(world)) has no pixel coordinate"
    coords = _frame_to_transform(last(wcs.steps).frame)(world)
    for index in (length(wcs.steps) - 1):-1:1
        inverse = _inverse_transform_checked(wcs.steps[index].transform)
        # The numerical inverse only proves undefinedness through analytically invertible steps.
        inverse isa NoInverse && break
        next = inverse(coords)
        if any(isnan, next)
            description = "the inverse of the step from frame \"$(wcs.steps[index].frame.name)\" to \"$(wcs.steps[index + 1].frame.name)\""
            throw(DomainError(Tuple(world), _undefined_message(subject, description, inverse, coords)))
        end
        coords = next
    end
    throw(DomainError(Tuple(world), subject))
end
