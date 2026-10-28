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

# `origin=1` numbers the first pixel center 1 (FITS, Julia arrays), like the `origin` argument of astropy.wcs: every pixel
# axis shifts and world coordinates never do. A Bool offset adds a strong zero for `origin=0`, keeping -0.0 and the input type.
@inline _origin_offset(origin) = origin == 0 ? false : origin == 1 ? true : _throw_origin(origin)
@noinline _throw_origin(origin) = throw(ArgumentError("origin must be 0 or 1, got $origin"))

@inline function _in_domain(wcs::WCS{N}, pixels::StaticVector{N}; with_bounding_box=true, origin=0) where {N}
    pixels = pixels .- _origin_offset(origin)
    all(isfinite, pixels) || return false
    return !with_bounding_box || wcs.bounding_box === nothing || pixels in wcs.bounding_box
end

"""
    in_domain(wcs, pixels...; origin=0) -> Bool
    in_domain(wcs, pixels::AbstractMatrix; origin=0) -> BitVector

Test whether finite pixel coordinates belong to the WCS input domain. Matrix columns are coordinates. Pixels are
zero-based; `origin=1` takes them one-based instead, on every input axis.
"""
function in_domain(wcs::WCS{N}, coords::AbstractVector; origin=0) where {N}
    length(coords) == N || throw(DimensionMismatch("pixel has length $(length(coords)), expected $N"))
    T = _coordinate_float_type(eltype(coords))
    return _in_domain(wcs, SVector{N,T}(coords); origin)
end

in_domain(wcs::WCS, coords::Tuple{Vararg{Real}}; origin=0) = _in_domain(wcs, _coordinate_vector(coords); origin)
in_domain(wcs::WCS, coords::Real...; origin=0) = in_domain(wcs, coords; origin)
in_domain(wcs::WCS{0}; origin=0) = _in_domain(wcs, SVector{0,Float64}(); origin)

function in_domain(wcs::WCS{N}, pixels::AbstractMatrix; origin=0) where {N}
    size(pixels, 1) == N || throw(DimensionMismatch("pixel matrix has $(size(pixels, 1)) rows, expected $N"))
    T = _coordinate_float_type(eltype(pixels))
    result = BitVector(undef, size(pixels, 2))
    @inbounds for (index, column) in enumerate(axes(pixels, 2))
        result[index] = _in_domain(wcs, SVector{N,T}(view(pixels, :, column)); origin)
    end
    return result
end

# `strict=false` returns `fill_value` for coordinates without a result (outside the bounding box or
# undefined under the transforms); `strict=true` throws for them instead. NaN input propagates as NaN.
@inline function _pixel_to_world(wcs::WCS{Nin,Nout}, coords::StaticVector{Nin}; with_bounding_box=true, fill_value=NaN, strict=false, origin=0) where {Nin,Nout}
    T = _coordinate_float_type(eltype(coords))
    given = SVector{Nin,T}(coords)
    pixels = given .- _origin_offset(origin)
    finite = all(isfinite, pixels)
    if with_bounding_box && wcs.bounding_box !== nothing && finite && !(pixels in wcs.bounding_box)
        strict && _throw_outside_box(given)
        return _filled_coordinate(Val(Nout), T, fill_value)
    end
    world = wcs.forward(pixels)
    # Transforms return NaN where a coordinate has no image; finite input that yields NaN is undefined.
    any(isnan, world) || return world
    strict || return _filled_coordinate(Val(Nout), T, fill_value)
    finite || return world
    _throw_undefined(wcs, pixels, given)
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

_angle_radians(unit) = unit == "deg" ? deg2rad(1.0) : unit == "rad" ? 1.0 : throw(ArgumentError("pixel_area requires celestial axes in deg or rad, got $(repr(unit))"))

# The latitude's transform position and the radian conversions of a 2D celestial output frame, as
# (latitude, lon * lat scale, lat scale).
function _celestial_area_scales(wcs::WCS)
    frame = output_frame(wcs)
    pixel_n_dim(wcs) == 2 && frame isa CelestialFrame{2} || throw(ArgumentError("pixel_area requires 2 pixel axes and a 2D CelestialFrame output frame, got $(pixel_n_dim(wcs)) pixel axes and a $(nameof(typeof(frame))){$(frame_dimension(frame))}"))
    latitude = _axis_rank(frame.axes_order, 2)
    radians = map(_angle_radians, axis_units(frame))
    return latitude, radians[1] * radians[2], radians[latitude]
end

# The world coordinate and Jacobian at a pixel, or `nothing` where there is no result: outside the bounding box, or an
# undefined coordinate or Jacobian. strict=true throws instead, and NaN input propagates, as in pixel_to_world.
# `transforms` evaluates part of the forward pipeline instead, such as the steps before a sky rotation.
@inline function _world_and_jacobian(wcs::WCS{Nin}, given::SVector{Nin}, transforms=wcs.forward; with_bounding_box=true, strict=false, origin=0) where {Nin}
    pixels = given .- _origin_offset(origin)
    finite = all(isfinite, pixels)
    if with_bounding_box && wcs.bounding_box !== nothing && finite && !(pixels in wcs.bounding_box)
        strict && _throw_outside_box(given)
        return nothing
    end
    world, jacobian = _value_and_jacobian(transforms, pixels)
    any(isnan, world) || any(isnan, jacobian) || return world, jacobian
    strict || return nothing
    finite || return world, jacobian
    any(isnan, world) && _throw_undefined(wcs, pixels, given)
    throw(DomainError(Tuple(given), "the Jacobian at pixel $(Tuple(given)) is undefined"))
end

_filled_matrix(::Val{M}, ::Val{N}, ::Type{T}, fill_value) where {M,N,T} = SMatrix{M,N,T}(ntuple(_ -> T(fill_value), Val(M * N)))

@inline function _pixel_to_world_jacobian(wcs::WCS{Nin,Nout}, coords::StaticVector{Nin}; fill_value=NaN, kwargs...) where {Nin,Nout}
    T = _coordinate_float_type(eltype(coords))
    result = _world_and_jacobian(wcs, SVector{Nin,T}(coords); kwargs...)
    return result === nothing ? _filled_matrix(Val(Nout), Val(Nin), T, fill_value) : last(result)
end

"""
    pixel_to_world_jacobian(wcs, pixels...; with_bounding_box=true, fill_value=NaN, strict=false, origin=0) -> SMatrix
    pixel_to_world_jacobian(wcs, pixels::AbstractMatrix; kwargs...) -> Vector{<:SMatrix}

Return the Jacobian of `pixel_to_world` at a pixel, the `world_n_dim × pixel_n_dim` matrix of `∂world_i/∂pixel_j`.
World coordinates are in transform order and in their frame's units, so a celestial WCS gives degrees per pixel with no
cos(latitude) factor. The coordinate forms and keywords act as in `pixel_to_world`: a pixel without a result, including
one whose Jacobian is undefined, gives a matrix of `fill_value`, or a `DomainError` with `strict=true`. Matrix columns
are pixels.
"""
function pixel_to_world_jacobian(wcs::WCS{Nin}, coords::AbstractVector; kwargs...) where {Nin}
    length(coords) == Nin || throw(DimensionMismatch("pixel has length $(length(coords)), expected $Nin"))
    return _pixel_to_world_jacobian(wcs, SVector{Nin,_coordinate_float_type(eltype(coords))}(coords); kwargs...)
end

pixel_to_world_jacobian(wcs::WCS, coords::Tuple{Vararg{Real}}; kwargs...) = _pixel_to_world_jacobian(wcs, _coordinate_vector(coords); kwargs...)
pixel_to_world_jacobian(wcs::WCS, coords::Real...; kwargs...) = pixel_to_world_jacobian(wcs, coords; kwargs...)

function pixel_to_world_jacobian(wcs::WCS{Nin}, pixels::AbstractMatrix; kwargs...) where {Nin}
    size(pixels, 1) == Nin || throw(DimensionMismatch("pixel matrix has $(size(pixels, 1)) rows, expected $Nin"))
    T = _coordinate_float_type(eltype(pixels))
    return map(column -> _pixel_to_world_jacobian(wcs, SVector{Nin,T}(view(pixels, :, column)); kwargs...), axes(pixels, 2))
end

# The Jacobian of the map world_to_pixel computes: the analytical inverse's, or the inverse of the forward Jacobian at
# the numerical inverse's pixel. Undefined results are thrown or filled before the bounding box applies, as there.
@inline function _world_to_pixel_jacobian(wcs::WCS{Nin,Nout}, coords::StaticVector{Nout}; with_bounding_box=true, fill_value=NaN, strict=false, origin=0, kwargs...) where {Nin,Nout}
    T = _coordinate_float_type(eltype(coords))
    world = SVector{Nout,T}(coords)
    offset = _origin_offset(origin)
    pixels, jacobian = if wcs.backward isa NoInverse
        solution = _numerical_solve(wcs, world; offset, kwargs...)
        solution, inv(_transform_jacobian_checked(wcs.forward, solution))
    else
        _value_and_jacobian(wcs.backward, world)
    end
    filled = _filled_matrix(Val(Nin), Val(Nout), T, fill_value)
    if any(isnan, pixels)
        strict || return filled
        all(isfinite, world) && _throw_undefined_inverse(wcs, world)
        return jacobian
    end
    if with_bounding_box && wcs.bounding_box !== nothing && !(pixels in wcs.bounding_box)
        strict && _throw_outside_box_inverse(world, pixels .+ offset)
        return filled
    end
    any(isnan, jacobian) || return jacobian
    strict || return filled
    throw(DomainError(Tuple(world), "the inverse Jacobian at world coordinate $(Tuple(world)) is undefined"))
end

"""
    world_to_pixel_jacobian(wcs, world...; with_bounding_box=true, fill_value=NaN, strict=false, origin=0, kwargs...) -> SMatrix
    world_to_pixel_jacobian(wcs, world::AbstractMatrix; kwargs...) -> Vector{<:SMatrix}

Return the Jacobian of `world_to_pixel` at a world coordinate, the `pixel_n_dim × world_n_dim` matrix of
`∂pixel_i/∂world_j`, with the coordinate forms, keywords, and units of `world_to_pixel`. With an analytical inverse it is
that inverse's Jacobian. Without one it is the inverse of the forward Jacobian at the pixel the numerical inverse finds,
which takes `world_to_pixel`'s solver keywords.

An analytical inverse that approximates the forward transform, such as a fitted polynomial, has a Jacobian that differs
from the inverse of the forward Jacobian by that approximation. The Jacobian of the exact inverse of the forward
transform is `inv(pixel_to_world_jacobian(wcs, numerical_inverse(wcs, world)))`.
"""
function world_to_pixel_jacobian(wcs::WCS{Nin,Nout}, coords::AbstractVector; kwargs...) where {Nin,Nout}
    length(coords) == Nout || throw(DimensionMismatch("world coordinate has length $(length(coords)), expected $Nout"))
    return _world_to_pixel_jacobian(wcs, SVector{Nout,_coordinate_float_type(eltype(coords))}(coords); kwargs...)
end

world_to_pixel_jacobian(wcs::WCS, coords::Tuple{Vararg{Real}}; kwargs...) = _world_to_pixel_jacobian(wcs, _coordinate_vector(coords); kwargs...)
world_to_pixel_jacobian(wcs::WCS, coords::Real...; kwargs...) = world_to_pixel_jacobian(wcs, coords; kwargs...)

function world_to_pixel_jacobian(wcs::WCS{Nin,Nout}, worlds::AbstractMatrix; kwargs...) where {Nin,Nout}
    size(worlds, 1) == Nout || throw(DimensionMismatch("world matrix has $(size(worlds, 1)) rows, expected $Nout"))
    T = _coordinate_float_type(eltype(worlds))
    return map(column -> _world_to_pixel_jacobian(wcs, SVector{Nout,T}(view(worlds, :, column)); kwargs...), axes(worlds, 2))
end

# A pipeline that ends in a rotation of the sphere, as Roman's and JWST's v2v3-to-world steps do, covers the same solid
# angle before it, where the longitude and latitude in degrees take fewer steps to differentiate and have not yet been
# rounded by the rotation. Returns the steps before it and the rotation, or `nothing`.
@generated function _sphere_rotation_split(transforms::Tuple)
    types = transforms.parameters
    n = length(types)
    n > 3 && types[n - 2] <: SphericalToCartesian && types[n - 1] <: Rotation3D && types[n] <: CartesianToSpherical || return :(nothing)
    Nin = supertype(types[1]).parameters[1]
    prefix = Expr(:tuple, (:(transforms[$i]) for i in 1:(n - 3))...)
    return :((ComposedTransform{$Nin,2,$(Tuple{types[1:(n - 3)]...})}($prefix), transforms[$(n - 1)]))
end

_sphere_rotation_split(transform::ComposedTransform) = _sphere_rotation_split(transform.transforms)
_sphere_rotation_split(::AbstractTransform) = nothing
# A Rotation3D built from a matrix need not be orthogonal, and only an orthogonal one preserves area.
_is_rotation(rotation::Rotation3D) = maximum(abs, rotation.matrix' * rotation.matrix - I) <= 1e-12

@inline function _pixel_area(wcs::WCS{2,2}, coords::StaticVector{2}, scales; fill_value=NaN, strict=false, kwargs...)
    T = _coordinate_float_type(eltype(coords))
    given = SVector{2,T}(coords)
    split = _sphere_rotation_split(wcs.forward)
    if split !== nothing && _is_rotation(last(split))
        result = _world_and_jacobian(wcs, given, first(split); strict, kwargs...)
        result === nothing && return T(fill_value)
        before, jacobian = result
        area = abs(det(jacobian)) * deg2rad(one(eltype(before)))^2 * cos(deg2rad(before[2]))
    else
        result = _world_and_jacobian(wcs, given; strict, kwargs...)
        result === nothing && return T(fill_value)
        world, jacobian = result
        latitude, solid_angle, latitude_radians = scales
        area = abs(det(jacobian)) * solid_angle * cos(world[latitude] * latitude_radians)
    end
    # A defined Jacobian can still give no finite area, as at a celestial pole, where the longitude derivative is infinite.
    isfinite(area) && return area
    strict || return T(fill_value)
    all(isfinite, given) || return area
    throw(DomainError(Tuple(given), "the pixel area at pixel $(Tuple(given)) is undefined"))
end

_unit_vector(lon, lat) = (sincos_lon = sincos(lon); sincos_lat = sincos(lat); SVector(sincos_lat[2] * sincos_lon[2], sincos_lat[2] * sincos_lon[1], sincos_lat[1]))
_cross2(u, v) = u[1] * v[2] - u[2] * v[1]
# The solid angle of a pixel from its corners a, b, c, d in order, from the cross product of its diagonals: of unit
# vectors projected on the pixel's direction, or of longitude and latitude in degrees with the metric at the center.
_sphere_corner_area(a, b, c, d) = (direction = a + b + c + d; abs(dot(direction, cross(c - a, d - b))) / (2 * norm(direction)))
_degree_corner_area(a, b, c, d) = abs(_cross2(c - a, d - b)) / 2 * cos(deg2rad((a[2] + b[2] + c[2] + d[2]) / 4)) * deg2rad(1.0)^2

# Corner rows are evaluated once and shared by the pixels on either side, so a grid costs about one evaluation per pixel.
function _pixel_area_grid(wcs::WCS{2,2}, xs, ys, corner, corner_area; with_bounding_box=true, fill_value=NaN, strict=false, origin=0)
    offset = _origin_offset(origin)
    corner_xs = (first(xs) - offset - 0.5):(last(xs) - offset + 0.5)
    corner_row!(row, y) = (for (k, x) in enumerate(corner_xs); row[k] = corner(SVector(x, y)); end; row)
    below = corner_row!(Vector{typeof(corner(SVector(0.0, 0.0)))}(undef, length(corner_xs)), first(ys) - offset - 0.5)
    above = similar(below)
    areas = Matrix{Float64}(undef, length(xs), length(ys))
    for (j, y) in enumerate(ys)
        corner_row!(above, y - offset + 0.5)
        for (i, x) in enumerate(xs)
            pixels = SVector{2,Float64}(x - offset, y - offset)
            if with_bounding_box && wcs.bounding_box !== nothing && !(pixels in wcs.bounding_box)
                strict && _throw_outside_box(SVector{2,Float64}(x, y))
                areas[i, j] = fill_value
                continue
            end
            area = corner_area(below[i], below[i + 1], above[i + 1], above[i])
            if !isfinite(area)
                strict && (any(isnan, wcs.forward(pixels)) ? _throw_undefined(wcs, pixels, SVector{2,Float64}(x, y)) : throw(DomainError((x, y), "the pixel area at pixel $((x, y)) is undefined")))
                area = fill_value
            end
            areas[i, j] = area
        end
        below, above = above, below
    end
    return areas
end

"""
    pixel_area(wcs, pixels...; with_bounding_box=true, fill_value=NaN, strict=false, origin=0) -> Real
    pixel_area(wcs, pixels::AbstractMatrix; kwargs...) -> Vector
    pixel_area(wcs, xs::AbstractUnitRange, ys::AbstractUnitRange; kwargs...) -> Matrix

Return the solid angle in steradians that a pixel covers on the sky, `|det J| cos(lat)`, with `J` the Jacobian of
the world coordinates with respect to the pixel at its center. The WCS needs 2 pixel axes and a 2D `CelestialFrame`
output frame in `deg` or `rad`. Matrix columns are pixels, and the keywords act as in `pixel_to_world`.

Given ranges of integer pixel coordinates, it returns the areas of the pixels in the grid `xs × ys` as a matrix
indexed like `data[x, y]`, so `pixel_area(wcs, axes(data)...; origin=1)` covers an image. Each is the solid angle within
the pixel's four corners, which neighboring pixels share, so the grid costs about one evaluation per pixel. Integrated
over the pixel rather than taken at its center, it differs from the per-pixel area by the variation of the area across
the pixel, about 1e-10 for Roman.
"""
function pixel_area(wcs::WCS, coords::AbstractVector; kwargs...)
    scales = _celestial_area_scales(wcs)
    length(coords) == 2 || throw(DimensionMismatch("pixel has length $(length(coords)), expected 2"))
    return _pixel_area(wcs, SVector{2,_coordinate_float_type(eltype(coords))}(coords), scales; kwargs...)
end

pixel_area(wcs::WCS, coords::Tuple{Vararg{Real}}; kwargs...) = pixel_area(wcs, _coordinate_vector(coords); kwargs...)
pixel_area(wcs::WCS, coords::Real...; kwargs...) = pixel_area(wcs, coords; kwargs...)

function pixel_area(wcs::WCS, xs::AbstractUnitRange{<:Integer}, ys::AbstractUnitRange{<:Integer}; kwargs...)
    latitude, _, latitude_radians = _celestial_area_scales(wcs)
    split = _sphere_rotation_split(wcs.forward)
    split !== nothing && _is_rotation(last(split)) && return _pixel_area_grid(wcs, xs, ys, first(split), _degree_corner_area; kwargs...)
    longitude_radians = _angle_radians(axis_units(output_frame(wcs))[3 - latitude])
    corner(pixels) = (world = wcs.forward(pixels); _unit_vector(world[3 - latitude] * longitude_radians, world[latitude] * latitude_radians))
    return _pixel_area_grid(wcs, xs, ys, corner, _sphere_corner_area; kwargs...)
end

function pixel_area(wcs::WCS, pixels::AbstractMatrix; kwargs...)
    scales = _celestial_area_scales(wcs)
    size(pixels, 1) == 2 || throw(DimensionMismatch("pixel matrix has $(size(pixels, 1)) rows, expected 2"))
    T = _coordinate_float_type(eltype(pixels))
    return map(column -> _pixel_area(wcs, SVector{2,T}(view(pixels, :, column)), scales; kwargs...), axes(pixels, 2))
end

@inline function _world_to_pixel(wcs::WCS{Nin,Nout}, coords::StaticVector{Nout}; with_bounding_box=true, fill_value=NaN, strict=false, origin=0, kwargs...) where {Nin,Nout}
    T = _coordinate_float_type(eltype(coords))
    world = SVector{Nout,T}(coords)
    offset = _origin_offset(origin)
    pixels = wcs.backward isa NoInverse ? _numerical_solve(wcs, world; offset, kwargs...) : wcs.backward(world)
    return _finish_inverse(wcs, world, pixels; with_bounding_box, fill_value, strict, offset)
end

# Shared by the analytical and numerical inverses: undefined results are thrown or filled before the bounding box applies.
# `pixels` are zero-based; `offset` shifts only a returned pixel, never `fill_value`.
@inline function _finish_inverse(wcs::WCS{Nin}, world, pixels; with_bounding_box, fill_value, strict, offset=false) where {Nin}
    if any(isnan, pixels)
        strict || return _filled_coordinate(Val(Nin), eltype(pixels), fill_value)
        all(isfinite, world) && _throw_undefined_inverse(wcs, world)
        return pixels
    end
    if with_bounding_box && wcs.bounding_box !== nothing && !(pixels in wcs.bounding_box)
        strict && _throw_outside_box_inverse(world, pixels .+ offset)
        return _filled_coordinate(Val(Nin), eltype(pixels), fill_value)
    end
    return pixels .+ offset
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

@inline function _in_image(wcs::WCS{Nin,Nout}, world::StaticVector{Nout}; with_bounding_box=true, tolerance=1e-5, maxiter=30, initial_guess=nothing, max_backtracks=12, origin=0) where {Nin,Nout}
    target = SVector{Nout,_coordinate_float_type(eltype(world))}(world)
    converged, pixels = try
        if wcs.backward isa NoInverse
            # The convergence flag replaces a thrown NoConvergence, which would cost far more than the solve.
            converged, pixels, _, _ = _numerical_result(wcs, target; tolerance, maxiter, initial_guess, max_backtracks, offset=_origin_offset(origin))
            (converged, pixels)
        else
            (true, wcs.backward(target))
        end
    catch error
        # DomainError remains only for zero-dimensional inverses, which cannot carry NaN.
        error isa DomainError || rethrow()
        return false
    end
    return converged && _in_domain(wcs, pixels; with_bounding_box)
end

"""
    in_image(wcs, world...; with_bounding_box=true, kwargs...) -> Bool
    in_image(wcs, worlds::AbstractMatrix; with_bounding_box=true, kwargs...) -> BitVector

Test whether world coordinates invert successfully to finite pixels in the WCS input domain. Matrix columns are coordinates.
`origin=1` takes `initial_guess` as one-based pixels.
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

# `given` is the caller's pixel, which differs from the zero-based `pixels` when `origin=1`.
@noinline function _throw_undefined(wcs::WCS, pixels, given=pixels)
    subject = "pixel $(Tuple(given)) has no world coordinate"
    coords = pixels
    for index in 1:(length(wcs.steps) - 1)
        step = wcs.steps[index]
        next = step.transform(coords)
        if any(isnan, next)
            description = "the step from frame \"$(step.frame.name)\" to \"$(wcs.steps[index + 1].frame.name)\""
            throw(DomainError(Tuple(given), _undefined_message(subject, description, step.transform, coords)))
        end
        coords = next
    end
    throw(DomainError(Tuple(given), subject))
end

@noinline function _throw_undefined_inverse(wcs::WCS, world)
    subject = "world coordinate $(Tuple(world)) has no pixel coordinate"
    coords = world
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
