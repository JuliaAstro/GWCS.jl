"""Construction-time description of a spherical projection algorithm."""
abstract type AbstractSphericalProjection end

"""Construction-time description of a zenithal spherical projection."""
abstract type AbstractZenithalProjection <: AbstractSphericalProjection end

"""The parameter-free TAN/gnomonic projection."""
struct GnomonicProjection <: AbstractZenithalProjection end

"""The parameter-free STG/stereographic projection."""
struct StereographicProjection <: AbstractZenithalProjection end

"""The parameter-free ARC/zenithal-equidistant projection."""
struct ZenithalEquidistant <: AbstractZenithalProjection end

"""The parameter-free ZEA/zenithal-equal-area projection."""
struct ZenithalEqualArea <: AbstractZenithalProjection end

"""The AZP/zenithal-perspective projection with distance `mu` and tilt `gamma`."""
struct ZenithalPerspective{T<:AbstractFloat} <: AbstractZenithalProjection
    mu::T
    gamma::T
    scale::T
    tan_gamma::T
    sec_gamma::T
    cos_gamma::T
    sin_gamma::T
    divergence_latitude::T
    mu_cos_gamma::T
    function ZenithalPerspective{T}(mu, gamma, scale, tan_gamma, sec_gamma, cos_gamma, sin_gamma, divergence_latitude, mu_cos_gamma) where {T}
        return new{T}(mu, gamma, scale, tan_gamma, sec_gamma, cos_gamma, sin_gamma, divergence_latitude, mu_cos_gamma)
    end
end

function ZenithalPerspective{T}(mu, gamma) where {T}
    mu_value = T(mu)
    gamma_value = T(gamma)
    all(isfinite, (mu_value, gamma_value)) || throw(ArgumentError("zenithal perspective parameters must be finite"))
    scale = rad2deg(one(T)) * (mu_value + one(T))
    iszero(scale) && throw(ArgumentError("zenithal perspective projection is undefined for mu = -1"))
    sin_gamma = sind(gamma_value)
    cos_gamma = cosd(gamma_value)
    iszero(cos_gamma) && throw(ArgumentError("zenithal perspective projection is undefined when cos(gamma) = 0"))
    tan_gamma = sin_gamma / cos_gamma
    sec_gamma = inv(cos_gamma)
    divergence_latitude = abs(mu_value) > one(T) ? rad2deg(asin(-inv(mu_value))) : T(-90)
    return ZenithalPerspective{T}(mu_value, gamma_value, scale, tan_gamma, sec_gamma, cos_gamma, sin_gamma, divergence_latitude, mu_value * cos_gamma)
end

function ZenithalPerspective(mu::Real=0, gamma::Real=0)
    T = promote_type(typeof(float(mu)), typeof(float(gamma)))
    return ZenithalPerspective{T}(mu, gamma)
end

"""The SZP/slant-zenithal-perspective projection."""
struct SlantZenithalPerspective{T<:AbstractFloat} <: AbstractZenithalProjection
    mu::T
    phi0::T
    theta0::T
    xp::T
    yp::T
    zp::T
    divergence_latitude::T
    function SlantZenithalPerspective{T}(mu, phi0, theta0, xp, yp, zp, divergence_latitude) where {T}
        return new{T}(mu, phi0, theta0, xp, yp, zp, divergence_latitude)
    end
end

function SlantZenithalPerspective{T}(mu, phi0, theta0) where {T}
    mu_value = T(mu)
    phi0_value = T(phi0)
    theta0_value = T(theta0)
    all(isfinite, (mu_value, phi0_value, theta0_value)) || throw(ArgumentError("slant zenithal perspective parameters must be finite"))
    -T(90) <= theta0_value <= T(90) || throw(ArgumentError("slant zenithal perspective theta0 must be between -90 and 90 degrees"))
    sin_phi0, cos_phi0 = sind(phi0_value), cosd(phi0_value)
    sin_theta0, cos_theta0 = sind(theta0_value), cosd(theta0_value)
    xp = -mu_value * cos_theta0 * sin_phi0
    yp = mu_value * cos_theta0 * cos_phi0
    zp = mu_value * sin_theta0 + one(T)
    iszero(zp) && throw(ArgumentError("slant zenithal perspective projection requires mu * sin(theta0) + 1 != 0"))
    divergence_latitude = abs(zp - one(T)) < one(T) ? rad2deg(asin(one(T) - zp)) : T(-90)
    return SlantZenithalPerspective{T}(mu_value, phi0_value, theta0_value, xp, yp, zp, divergence_latitude)
end

function SlantZenithalPerspective(mu::Real=0, phi0::Real=0, theta0::Real=90)
    T = promote_type(typeof(float(mu)), typeof(float(phi0)), typeof(float(theta0)))
    return SlantZenithalPerspective{T}(mu, phi0, theta0)
end

"""The SIN/slant-orthographic projection with obliqueness parameters `xi` and `eta`."""
struct SlantOrthographic{T<:AbstractFloat} <: AbstractZenithalProjection
    xi::T
    eta::T
    squared_slant::T
    function SlantOrthographic{T}(xi, eta, squared_slant) where {T}
        return new{T}(xi, eta, squared_slant)
    end
end

function SlantOrthographic{T}(xi, eta) where {T}
    xi_value = T(xi)
    eta_value = T(eta)
    all(isfinite, (xi_value, eta_value)) || throw(ArgumentError("slant orthographic parameters must be finite"))
    return SlantOrthographic{T}(xi_value, eta_value, xi_value^2 + eta_value^2)
end

function SlantOrthographic(xi::Real=0, eta::Real=0)
    T = promote_type(typeof(float(xi)), typeof(float(eta)))
    return SlantOrthographic{T}(xi, eta)
end

"""The AIR/Airy projection with error-minimization latitude `theta_b`."""
struct AiryProjection{T<:AbstractFloat} <: AbstractZenithalProjection
    theta_b::T
    coefficient::T
    slope::T
    function AiryProjection{T}(theta_b, coefficient, slope) where {T}
        return new{T}(theta_b, coefficient, slope)
    end
end

function AiryProjection{T}(theta_b) where {T}
    theta = T(theta_b)
    isfinite(theta) || throw(ArgumentError("Airy theta_b must be finite"))
    theta > -oftype(theta, 90) || throw(ArgumentError("Airy theta_b must be greater than -90 degrees"))
    xi_b = deg2rad((oftype(theta, 90) - theta) / 2)
    coefficient = iszero(xi_b) ? -oftype(theta, 0.5) : log(cos(xi_b)) / tan(xi_b)^2
    slope = oftype(theta, 0.5) - coefficient
    return AiryProjection{typeof(theta)}(theta, coefficient, slope)
end

AiryProjection(theta_b::Real=90) = AiryProjection{typeof(float(theta_b))}(theta_b)

"""Convert projection-plane coordinates in degrees to native spherical coordinates in degrees."""
struct IntermediateToNative{P<:AbstractSphericalProjection} <: AbstractTransform{2,2}
    projection::P
end

"""Convert native spherical coordinates in degrees to projection-plane coordinates in degrees."""
struct NativeToIntermediate{P<:AbstractSphericalProjection} <: AbstractTransform{2,2}
    projection::P
end

IntermediateToNative{P}() where {P<:AbstractSphericalProjection} = IntermediateToNative(P())
NativeToIntermediate{P}() where {P<:AbstractSphericalProjection} = NativeToIntermediate(P())

inverse_transform(transform::IntermediateToNative) = NativeToIntermediate(transform.projection)
inverse_transform(transform::NativeToIntermediate) = IntermediateToNative(transform.projection)

@inline function (transform::IntermediateToNative)(coords::StaticVector{2})
    T = _evaluation_type(eltype(coords))
    phi, theta = _intermediate_to_native(transform.projection, T(coords[1]), T(coords[2]))
    return SVector(rad2deg(phi), rad2deg(theta))
end

@inline function (transform::NativeToIntermediate)(coords::StaticVector{2})
    T = _evaluation_type(eltype(coords))
    phi = deg2rad(T(coords[1]))
    theta = deg2rad(T(coords[2]))
    x, y = _native_to_intermediate(transform.projection, phi, theta)
    return SVector(x, y)
end

@inline function _zenithal_phi(x, y)
    iszero(x) && iszero(y) && return zero(promote_type(typeof(x), typeof(y)))
    return atan(x, -y)
end

@inline function _zenithal_xy(radius, phi)
    sine, cosine = sincos(phi)
    return radius * sine, -radius * cosine
end

@inline function _projection_tolerance(::Type{T}) where {T}
    return T(100) * eps(T)
end

function _native_latitude(theta, projection)
    limit = oftype(theta, pi / 2)
    tolerance = _projection_tolerance(typeof(theta))
    abs(theta) <= limit + tolerance || throw(DomainError(rad2deg(theta), "$(typeof(projection)) native latitude must be between -90 and 90 degrees"))
    return clamp(theta, -limit, limit)
end

@inline function _intermediate_to_native(::GnomonicProjection, x, y)
    T = promote_type(typeof(x), typeof(y))
    return _zenithal_phi(x, y), atan(rad2deg(one(T)), hypot(x, y))
end

@inline function _native_to_intermediate(projection::GnomonicProjection, phi, theta)
    theta = _native_latitude(theta, projection)
    sine, cosine = sincos(theta)
    sine > zero(sine) || throw(DomainError(rad2deg(theta), "TAN native latitude must be positive"))
    return _zenithal_xy(rad2deg(cosine / sine), phi)
end

@inline function _intermediate_to_native(::StereographicProjection, x, y)
    T = promote_type(typeof(x), typeof(y))
    radius = deg2rad(hypot(x, y)) / 2
    sine = (one(T) - radius^2) / (one(T) + radius^2)
    return _zenithal_phi(x, y), asin(clamp(sine, -one(T), one(T)))
end

@inline function _native_to_intermediate(projection::StereographicProjection, phi, theta)
    theta = _native_latitude(theta, projection)
    denominator = one(theta) + sin(theta)
    iszero(denominator) && throw(DomainError(rad2deg(theta), "STG is singular at native latitude -90 degrees"))
    return _zenithal_xy(rad2deg(2cos(theta) / denominator), phi)
end

@inline function _intermediate_to_native(::ZenithalEquidistant, x, y)
    T = promote_type(typeof(x), typeof(y))
    theta = oftype(one(T), pi / 2) - deg2rad(hypot(x, y))
    theta >= -oftype(theta, pi / 2) - _projection_tolerance(T) || throw(DomainError(SVector(x, y), "ARC intermediate coordinate is outside the projection domain"))
    return _zenithal_phi(x, y), max(theta, -oftype(theta, pi / 2))
end

@inline function _native_to_intermediate(projection::ZenithalEquidistant, phi, theta)
    theta = _native_latitude(theta, projection)
    return _zenithal_xy(rad2deg(oftype(theta, pi / 2) - theta), phi)
end

@inline function _intermediate_to_native(::ZenithalEqualArea, x, y)
    T = promote_type(typeof(x), typeof(y))
    argument = deg2rad(hypot(x, y)) / 2
    tolerance = _projection_tolerance(T)
    argument <= one(T) + tolerance || throw(DomainError(SVector(x, y), "ZEA intermediate coordinate is outside the projection domain"))
    theta = oftype(one(T), pi / 2) - 2asin(clamp(argument, zero(T), one(T)))
    return _zenithal_phi(x, y), theta
end

@inline function _native_to_intermediate(projection::ZenithalEqualArea, phi, theta)
    theta = _native_latitude(theta, projection)
    radius = rad2deg(2sin((oftype(theta, pi / 2) - theta) / 2))
    return _zenithal_xy(radius, phi)
end

function _intermediate_to_native(projection::SlantOrthographic, x, y)
    T = promote_type(typeof(x), typeof(y))
    xi = T(projection.xi)
    eta = T(projection.eta)
    squared_slant = T(projection.squared_slant)
    xr = deg2rad(x)
    yr = deg2rad(y)
    radius_squared = xr^2 + yr^2

    if iszero(squared_slant)
        radius_squared <= one(T) + _projection_tolerance(T) || throw(DomainError(SVector(x, y), "SIN intermediate coordinate is outside the projection domain"))
        radius_squared = min(radius_squared, one(T))
        phi = iszero(radius_squared) ? zero(T) : atan(xr, -yr)
        theta = radius_squared < T(0.5) ? acos(sqrt(radius_squared)) : asin(sqrt(one(T) - radius_squared))
        return phi, theta
    end

    xy = xr * xi + yr * eta
    if radius_squared < T(1e-10)
        z = radius_squared / 2
        theta = oftype(one(T), pi / 2) - sqrt(radius_squared / (one(T) + xy))
    else
        a = squared_slant + one(T)
        b = xy - squared_slant
        c = radius_squared - 2xy + squared_slant - one(T)
        discriminant = b^2 - a * c
        tolerance = _projection_tolerance(T)
        discriminant >= -tolerance || throw(DomainError(SVector(x, y), "SIN intermediate coordinate has no native solution"))
        root = sqrt(max(discriminant, zero(T)))
        sine1 = (-b + root) / a
        sine2 = (-b - root) / a
        sine = max(sine1, sine2)
        if sine > one(T) + tolerance
            sine = min(sine1, sine2)
        end
        -one(T) - tolerance <= sine <= one(T) + tolerance || throw(DomainError(SVector(x, y), "SIN intermediate coordinate has no native solution"))
        sine = clamp(sine, -one(T), one(T))
        theta = asin(sine)
        z = one(T) - sine
    end
    dx = xr - xi * z
    dy = yr - eta * z
    phi = iszero(dx) && iszero(dy) ? zero(T) : atan(dx, -dy)
    return phi, theta
end

function _native_to_intermediate(projection::SlantOrthographic, phi, theta)
    T = promote_type(typeof(phi), typeof(theta))
    theta = _native_latitude(theta, projection)
    xi = T(projection.xi)
    eta = T(projection.eta)
    sine_phi, cosine_phi = sincos(phi)
    sine_theta, cosine_theta = sincos(theta)
    if iszero(T(projection.squared_slant))
        theta >= zero(T) || throw(DomainError(rad2deg(theta), "standard SIN native latitude must be non-negative"))
    else
        minimum_theta = -atan(xi * sine_phi - eta * cosine_phi)
        theta >= minimum_theta - _projection_tolerance(T) || throw(DomainError(rad2deg(theta), "slant SIN native coordinate is outside the projection domain"))
    end
    offset = one(T) - sine_theta
    return rad2deg(cosine_theta * sine_phi + xi * offset), rad2deg(-cosine_theta * cosine_phi + eta * offset)
end

function _intermediate_to_native(projection::ZenithalPerspective, x, y)
    T = promote_type(typeof(x), typeof(y))
    mu = T(projection.mu)
    cosine_gamma = T(projection.cos_gamma)
    sine_gamma = T(projection.sin_gamma)
    yc = y * cosine_gamma
    radius = hypot(x, yc)
    iszero(radius) && return zero(T), oftype(one(T), pi / 2)
    q = T(projection.scale) + y * sine_gamma
    iszero(q) && throw(DomainError(SVector(x, y), "AZP intermediate coordinate is singular"))
    ratio = radius / q
    argument = ratio * mu / sqrt(ratio^2 + one(T))
    tolerance = _projection_tolerance(T)
    abs(argument) <= one(T) + tolerance || throw(DomainError(SVector(x, y), "AZP intermediate coordinate has no native solution"))
    angle = atan(one(T), ratio)
    offset = asin(clamp(argument, -one(T), one(T)))
    first_solution = angle - offset
    second_solution = angle + offset + oftype(angle, pi)
    first_solution > oftype(angle, pi / 2) && (first_solution -= oftype(angle, 2pi))
    second_solution > oftype(angle, pi / 2) && (second_solution -= oftype(angle, 2pi))
    theta = max(first_solution, second_solution)
    return atan(x, -yc), theta
end

function _native_to_intermediate(projection::ZenithalPerspective, phi, theta)
    T = promote_type(typeof(phi), typeof(theta))
    theta = _native_latitude(theta, projection)
    mu = T(projection.mu)
    sine_phi, cosine_phi = sincos(phi)
    sine_theta, cosine_theta = sincos(theta)
    tangent_gamma = T(projection.tan_gamma)
    denominator = mu + sine_theta + cosine_theta * tangent_gamma * cosine_phi
    iszero(denominator) && throw(DomainError(rad2deg(theta), "AZP native coordinate is singular"))

    theta_degrees = rad2deg(theta)
    theta_degrees >= T(projection.divergence_latitude) - _projection_tolerance(T) || throw(DomainError(theta_degrees, "AZP native coordinate overlaps another projection branch"))
    if abs(T(projection.mu_cos_gamma)) < one(T)
        slope = tangent_gamma * cosine_phi
        argument = mu / sqrt(one(T) + slope^2)
        if abs(argument) <= one(T)
            angle = atan(-slope)
            offset = asin(argument)
            first_limit = angle - offset
            second_limit = angle + offset + oftype(angle, pi)
            first_limit > oftype(angle, pi / 2) && (first_limit -= oftype(angle, 2pi))
            second_limit > oftype(angle, pi / 2) && (second_limit -= oftype(angle, 2pi))
            theta >= max(first_limit, second_limit) - _projection_tolerance(T) || throw(DomainError(theta_degrees, "AZP native coordinate is beyond the projection divergence"))
        end
    end

    radius = T(projection.scale) * cosine_theta / denominator
    return radius * sine_phi, -radius * cosine_phi * T(projection.sec_gamma)
end

function _intermediate_to_native(projection::SlantZenithalPerspective, x, y)
    T = promote_type(typeof(x), typeof(y))
    inverse_scale = deg2rad(one(T))
    xr = x * inverse_scale
    yr = y * inverse_scale
    xp = T(projection.xp)
    yp = T(projection.yp)
    zp = T(projection.zp)
    radius_squared = xr^2 + yr^2
    x1 = (xr - xp) / zp
    y1 = (yr - yp) / zp
    xy = xr * x1 + yr * y1

    if radius_squared < T(1e-10)
        z = radius_squared / 2
        theta = oftype(one(T), pi / 2) - sqrt(radius_squared / (one(T) + xy))
    else
        t = x1^2 + y1^2
        a = t + one(T)
        b = xy - t
        c = radius_squared - 2xy + t - one(T)
        discriminant = b^2 - a * c
        tolerance = _projection_tolerance(T)
        discriminant >= -tolerance || throw(DomainError(SVector(x, y), "SZP intermediate coordinate has no native solution"))
        root = sqrt(max(discriminant, zero(T)))
        sine1 = (-b + root) / a
        sine2 = (-b - root) / a
        sine = max(sine1, sine2)
        if sine > one(T) + tolerance
            sine = min(sine1, sine2)
        end
        -one(T) - tolerance <= sine <= one(T) + tolerance || throw(DomainError(SVector(x, y), "SZP intermediate coordinate has no native solution"))
        sine = clamp(sine, -one(T), one(T))
        theta = asin(sine)
        z = one(T) - sine
    end
    dx = xr - x1 * z
    dy = yr - y1 * z
    phi = iszero(dx) && iszero(dy) ? zero(T) : atan(dx, -dy)
    return phi, theta
end

function _native_to_intermediate(projection::SlantZenithalPerspective, phi, theta)
    T = promote_type(typeof(phi), typeof(theta))
    theta = _native_latitude(theta, projection)
    xp = T(projection.xp)
    yp = T(projection.yp)
    zp = T(projection.zp)
    sine_phi, cosine_phi = sincos(phi)
    sine_theta, cosine_theta = sincos(theta)
    offset = one(T) - sine_theta
    denominator = zp - offset
    iszero(denominator) && throw(DomainError(rad2deg(theta), "SZP native coordinate is singular"))

    theta_degrees = rad2deg(theta)
    theta_degrees >= T(projection.divergence_latitude) - _projection_tolerance(T) || throw(DomainError(theta_degrees, "SZP native coordinate is beyond the projection divergence"))
    mu = T(projection.mu)
    if abs(mu) > one(T)
        slope = xp * sine_phi - yp * cosine_phi
        radicand = (zp - one(T)) * zp - one(T) + slope^2
        if radicand > zero(T)
            argument = inv(sqrt(radicand))
            if abs(argument) <= one(T)
                angle = atan(slope, zp - one(T))
                offset_angle = asin(argument)
                first_limit = angle - offset_angle
                second_limit = angle + offset_angle + oftype(angle, pi)
                first_limit > oftype(angle, pi / 2) && (first_limit -= oftype(angle, 2pi))
                second_limit > oftype(angle, pi / 2) && (second_limit -= oftype(angle, 2pi))
                theta >= max(first_limit, second_limit) - _projection_tolerance(T) || throw(DomainError(theta_degrees, "SZP native coordinate overlaps another projection branch"))
            end
        end
    end

    scale = rad2deg(one(T))
    radius = scale * zp * cosine_theta / denominator
    u = scale * xp * offset / denominator
    v = scale * yp * offset / denominator
    return radius * sine_phi - u, -radius * cosine_phi - v
end

@inline function _airy_radius(projection::AiryProjection, theta)
    T = typeof(theta)
    xi = (oftype(theta, pi / 2) - theta) / 2
    xi >= zero(T) || throw(DomainError(rad2deg(theta), "Airy native latitude must not exceed 90 degrees"))
    xi < oftype(theta, pi / 2) || throw(DomainError(rad2deg(theta), "Airy projection is singular at native latitude -90 degrees"))
    coefficient = T(projection.coefficient)
    if xi < cbrt(eps(T))
        return rad2deg(2xi * T(projection.slope))
    end
    tangent = tan(xi)
    return -2rad2deg(one(T)) * (log(cos(xi)) / tangent + coefficient * tangent)
end

function _intermediate_to_native(projection::AiryProjection, x, y)
    T = promote_type(typeof(x), typeof(y))
    phi = _zenithal_phi(x, y)
    target = hypot(x, y)
    iszero(target) && return phi, oftype(one(T), pi / 2)

    lower = zero(T)
    upper = oftype(one(T), pi / 2) - sqrt(eps(T))
    upper_radius = _airy_radius(projection, oftype(one(T), pi / 2) - 2upper)
    target <= upper_radius || throw(DomainError(SVector(x, y), "Airy intermediate coordinate is outside the supported projection domain"))
    for _ in 1:96
        midpoint = (lower + upper) / 2
        theta = oftype(one(T), pi / 2) - 2midpoint
        radius = _airy_radius(projection, theta)
        if radius < target
            lower = midpoint
        else
            upper = midpoint
        end
        upper - lower <= T(8) * eps(T) * max(one(T), midpoint) && break
    end
    xi = (lower + upper) / 2
    return phi, oftype(one(T), pi / 2) - 2xi
end

function _native_to_intermediate(projection::AiryProjection, phi, theta)
    theta = _native_latitude(theta, projection)
    return _zenithal_xy(_airy_radius(projection, theta), phi)
end

function transform_jacobian(transform::IntermediateToNative{<:AiryProjection}, coords::StaticVector{2})
    native = transform(coords)
    inverse_jacobian = _forwarddiff_jacobian(NativeToIntermediate(transform.projection), native)
    return inv(inverse_jacobian)
end

"""
    native_fiducial(projection)

Return the native longitude and latitude in degrees of a projection's
fiducial point.
"""
native_fiducial(::AbstractZenithalProjection) = (0.0, 90.0)

"""Return the FITS Paper II default `LONPOLE` in degrees."""
function default_lonpole(reference_latitude::Real, native_latitude::Real, native_longitude::Real=0)
    values = promote(float(reference_latitude), float(native_latitude), float(native_longitude))
    all(isfinite, values) || throw(ArgumentError("celestial and native fiducial coordinates must be finite"))
    lonpole = values[1] < values[2] ? values[3] + oftype(values[3], 180) : values[3]
    wrapped = mod(lonpole, oftype(lonpole, 360))
    return wrapped > oftype(wrapped, 180) ? wrapped - oftype(wrapped, 360) : wrapped
end

function _native_pole(reference_longitude, reference_latitude, native_longitude, native_latitude, lonpole, latpole)
    values = promote(float(reference_longitude), float(reference_latitude), float(native_longitude), float(native_latitude), float(lonpole), float(latpole))
    all(isfinite, values) || throw(ArgumentError("celestial pole parameters must be finite"))
    T = typeof(values[1])
    longitude0, latitude0, phi0, theta0, phi_p, latitude_hint = values
    -T(90) <= latitude0 <= T(90) || throw(ArgumentError("reference latitude must be between -90 and 90 degrees"))
    -T(90) <= theta0 <= T(90) || throw(ArgumentError("native fiducial latitude must be between -90 and 90 degrees"))
    latitude_hint = clamp(latitude_hint, -T(90), T(90))
    tolerance = T(1e-10)

    if isapprox(theta0, T(90); atol=tolerance, rtol=zero(T))
        return longitude0, latitude0
    end

    sine_latitude0, cosine_latitude0 = sincos(deg2rad(latitude0))
    sine_theta0, cosine_theta0 = sincos(deg2rad(theta0))
    sine_phi, cosine_phi = sincos(deg2rad(phi_p - phi0))
    u = zero(T)
    v = zero(T)
    specified_by_latpole = false

    if iszero(phi_p - phi0)
        u = theta0
        v = T(90) - latitude0
    else
        x = cosine_theta0 * cosine_phi
        y = sine_theta0
        z = hypot(x, y)
        if iszero(z)
            iszero(sine_latitude0) || throw(DomainError(latitude0, "reference latitude must be zero for this degenerate pole configuration"))
            specified_by_latpole = true
        else
            sine_latitude_z = sine_latitude0 / z
            abs(sine_latitude_z) <= one(T) + tolerance || throw(DomainError(latitude0, "no native pole exists for the requested reference coordinate"))
            u = rad2deg(atan(y, x))
            v = rad2deg(acos(clamp(sine_latitude_z, -one(T), one(T))))
        end
    end

    latitude_pole = latitude_hint
    if !specified_by_latpole
        candidate1 = mod(u + v + T(180), T(360)) - T(180)
        candidate2 = mod(u - v + T(180), T(360)) - T(180)
        valid1 = abs(candidate1) <= T(90) + tolerance
        valid2 = abs(candidate2) <= T(90) + tolerance
        valid1 || valid2 || throw(DomainError((candidate1, candidate2), "no valid native-pole latitude exists"))
        latitude_pole = if valid1 && valid2
            abs(latitude_hint - candidate1) < abs(latitude_hint - candidate2) ? candidate1 : candidate2
        elseif valid1
            candidate1
        else
            candidate2
        end
        latitude_pole = clamp(latitude_pole, -T(90), T(90))
    end

    denominator_scale = cosd(latitude_pole) * cosine_latitude0
    longitude_pole = if abs(denominator_scale) < tolerance
        if abs(cosine_latitude0) < tolerance
            longitude0
        elseif latitude_pole > zero(T)
            longitude0 + phi_p - phi0 - T(180)
        else
            longitude0 - phi_p + phi0
        end
    else
        x = (sine_theta0 - sind(latitude_pole) * sine_latitude0) / denominator_scale
        y = sine_phi * cosine_theta0 / cosine_latitude0
        iszero(x) && iszero(y) && throw(DomainError((x, y), "native-pole longitude is indeterminate"))
        longitude0 - rad2deg(atan(y, x))
    end

    if longitude0 >= zero(T)
        longitude_pole < zero(T) && (longitude_pole += T(360))
        longitude_pole > T(360) && (longitude_pole -= T(360))
    else
        longitude_pole > zero(T) && (longitude_pole -= T(360))
        longitude_pole < -T(360) && (longitude_pole += T(360))
    end
    return longitude_pole, latitude_pole
end

"""
    celestial_transform(projection, reference; lonpole=nothing, latpole=90)

Construct an intermediate-plane-to-celestial transform using a two-coordinate
celestial reference point and FITS Paper II native-pole conventions. All
coordinates and pole parameters are in degrees.
"""
function celestial_transform(projection::AbstractSphericalProjection, reference; lonpole=nothing, latpole=90)
    reference isa Tuple || reference isa AbstractVector || throw(ArgumentError("celestial reference coordinate must be a tuple or vector"))
    length(reference) == 2 || throw(DimensionMismatch("celestial reference coordinate has length $(length(reference)), expected 2"))
    all(value -> value isa Real, reference) || throw(ArgumentError("celestial reference coordinates must be real numbers"))
    lonpole === nothing || lonpole isa Real || throw(ArgumentError("lonpole must be a real number or nothing"))
    latpole isa Real || throw(ArgumentError("latpole must be a real number"))
    phi0, theta0 = native_fiducial(projection)
    reference_longitude, reference_latitude = reference
    phi_p = lonpole === nothing ? default_lonpole(reference_latitude, theta0, phi0) : lonpole
    alpha_p, delta_p = _native_pole(reference_longitude, reference_latitude, phi0, theta0, phi_p, latpole)
    return IntermediateToNative(projection) | NativeToCelestial(alpha_p, delta_p, phi_p)
end
