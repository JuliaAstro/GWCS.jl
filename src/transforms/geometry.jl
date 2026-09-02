"""Convert longitude and latitude in degrees to a Cartesian unit vector."""
struct SphericalToCartesian <: AbstractTransform{2,3}
    wrap_lon_at::Int
    function SphericalToCartesian(wrap_lon_at)
        wrap_lon_at in (180, 360) || throw(ArgumentError("wrap_lon_at must be 180 or 360"))
        return new(Int(wrap_lon_at))
    end
end

SphericalToCartesian() = SphericalToCartesian(360)

@inline function (transform::SphericalToCartesian)(coords::StaticVector{2})
    lon = deg2rad(coords[1])
    lat = deg2rad(coords[2])
    sinlon, coslon = sincos(lon)
    sinlat, coslat = sincos(lat)
    return SVector(coslat * coslon, coslat * sinlon, sinlat)
end

inverse_transform(transform::SphericalToCartesian) = CartesianToSpherical(transform.wrap_lon_at)

"""Convert a Cartesian vector to longitude and latitude in degrees."""
struct CartesianToSpherical <: AbstractTransform{3,2}
    wrap_lon_at::Int
    function CartesianToSpherical(wrap_lon_at)
        wrap_lon_at in (180, 360) || throw(ArgumentError("wrap_lon_at must be 180 or 360"))
        return new(Int(wrap_lon_at))
    end
end

CartesianToSpherical() = CartesianToSpherical(360)

@inline function (transform::CartesianToSpherical)(coords::StaticVector{3})
    x, y, z = coords
    h = hypot(x, y)
    lon = rad2deg(atan(y, x))
    iszero(h) && (lon = zero(lon))
    transform.wrap_lon_at == 360 && (lon = mod(lon, oftype(lon, 360)))
    return SVector(lon, rad2deg(atan(z, h)))
end

inverse_transform(transform::CartesianToSpherical) = SphericalToCartesian(transform.wrap_lon_at)

"""Rotate native spherical coordinates in degrees to celestial coordinates."""
struct NativeToCelestial{T<:Real} <: AbstractTransform{2,2}
    alpha_p::T
    delta_p::T
    phi_p::T
end

function NativeToCelestial(alpha_p::Real, delta_p::Real, phi_p::Real)
    values = promote(alpha_p, delta_p, phi_p)
    return NativeToCelestial{typeof(values[1])}(values...)
end

@inline function (transform::NativeToCelestial)(coords::StaticVector{2})
    T = _evaluation_type(eltype(coords))
    phi = deg2rad(T(coords[1]))
    theta = deg2rad(T(coords[2]))
    alpha_p = deg2rad(T(transform.alpha_p))
    delta_p = deg2rad(T(transform.delta_p))
    phi_p = deg2rad(T(transform.phi_p))
    dphi = phi - phi_p
    sin_theta, cos_theta = sincos(theta)
    sin_dphi, cos_dphi = sincos(dphi)
    sin_delta_p, cos_delta_p = sincos(delta_p)
    delta = asin(clamp(sin_theta * sin_delta_p + cos_theta * cos_delta_p * cos_dphi, -one(T), one(T)))
    alpha = alpha_p + atan(-cos_theta * sin_dphi, sin_theta * cos_delta_p - cos_theta * sin_delta_p * cos_dphi)
    return SVector(mod(rad2deg(alpha), T(360)), rad2deg(delta))
end

inverse_transform(transform::NativeToCelestial) = CelestialToNative(transform.alpha_p, transform.delta_p, transform.phi_p)

"""Rotate celestial coordinates in degrees to native spherical coordinates."""
struct CelestialToNative{T<:Real} <: AbstractTransform{2,2}
    alpha_p::T
    delta_p::T
    phi_p::T
end

function CelestialToNative(alpha_p::Real, delta_p::Real, phi_p::Real)
    values = promote(alpha_p, delta_p, phi_p)
    return CelestialToNative{typeof(values[1])}(values...)
end

@inline function (transform::CelestialToNative)(coords::StaticVector{2})
    T = _evaluation_type(eltype(coords))
    alpha = deg2rad(T(coords[1]))
    delta = deg2rad(T(coords[2]))
    alpha_p = deg2rad(T(transform.alpha_p))
    delta_p = deg2rad(T(transform.delta_p))
    phi_p = deg2rad(T(transform.phi_p))
    delta_alpha = alpha - alpha_p
    sin_delta, cos_delta = sincos(delta)
    sin_delta_alpha, cos_delta_alpha = sincos(delta_alpha)
    sin_delta_p, cos_delta_p = sincos(delta_p)
    x = sin_delta * cos_delta_p - cos_delta * sin_delta_p * cos_delta_alpha
    y = -cos_delta * sin_delta_alpha
    z = sin_delta * sin_delta_p + cos_delta * cos_delta_p * cos_delta_alpha
    theta = abs(z) > T(0.99) ? copysign(acos(clamp(hypot(x, y), zero(T), one(T))), z) : asin(clamp(z, -one(T), one(T)))
    phi = mod(phi_p + atan(y, x) + T(pi), T(2pi)) - T(pi)
    return SVector(mod(rad2deg(phi), T(360)), rad2deg(theta))
end

inverse_transform(transform::CelestialToNative) = NativeToCelestial(transform.alpha_p, transform.delta_p, transform.phi_p)

@inline function _axis_rotation(angle, axis::Char)
    sine, cosine = sincos(deg2rad(angle))
    zero_angle = zero(angle)
    one_angle = one(angle)
    if axis == 'x'
        return SMatrix{3,3}(one_angle, zero_angle, zero_angle, zero_angle, cosine, -sine, zero_angle, sine, cosine)
    elseif axis == 'y'
        return SMatrix{3,3}(cosine, zero_angle, sine, zero_angle, one_angle, zero_angle, -sine, zero_angle, cosine)
    elseif axis == 'z'
        return SMatrix{3,3}(cosine, -sine, zero_angle, sine, cosine, zero_angle, zero_angle, zero_angle, one_angle)
    end
    throw(ArgumentError("rotation axis must be x, y, or z"))
end

"""Three-dimensional rotation represented by a precomputed static matrix."""
struct Rotation3D{T,L} <: AbstractTransform{3,3}
    matrix::SMatrix{3,3,T,L}
end

function Rotation3D(matrix::AbstractMatrix)
    size(matrix) == (3, 3) || throw(DimensionMismatch("rotation matrix must be 3 x 3"))
    T = eltype(matrix)
    return Rotation3D{T,9}(SMatrix{3,3,T}(matrix))
end

function Rotation3D(angles::AbstractVector, axes_order::AbstractString)
    length(angles) == length(axes_order) || throw(ArgumentError("the number of angles must match axes_order"))
    T = float(eltype(angles))
    matrix = SMatrix{3,3,T}(I)
    for (angle, axis) in zip(angles, axes_order)
        matrix = _axis_rotation(T(angle), axis) * matrix
    end
    return Rotation3D(matrix)
end

@inline function (transform::Rotation3D)(coords::StaticVector{3})
    T = _evaluation_type(eltype(coords))
    return SMatrix{3,3,T}(transform.matrix) * SVector{3,T}(coords)
end

inverse_transform(transform::Rotation3D) = Rotation3D(transpose(transform.matrix))
