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

"""The CAR/plate carree projection."""
struct PlateCarreeProjection <: AbstractSphericalProjection end

"""The CEA/cylindrical equal-area projection with parameter `lambda`."""
struct CylindricalEqualArea{T<:AbstractFloat} <: AbstractSphericalProjection
    lambda::T
    function CylindricalEqualArea{T}(lambda) where {T}
        value = T(lambda)
        isfinite(value) && 0 < value <= 1 || throw(ArgumentError("cylindrical equal-area lambda must be in (0, 1]"))
        return new{T}(value)
    end
end

CylindricalEqualArea(lambda::Real=1) = CylindricalEqualArea{typeof(float(lambda))}(lambda)

"""The CYP/cylindrical perspective projection."""
struct CylindricalPerspective{T<:AbstractFloat} <: AbstractSphericalProjection
    mu::T
    lambda::T
    function CylindricalPerspective{T}(mu, lambda) where {T}
        mu_value, lambda_value = T(mu), T(lambda)
        all(isfinite, (mu_value, lambda_value)) || throw(ArgumentError("cylindrical perspective parameters must be finite"))
        iszero(lambda_value) && throw(ArgumentError("cylindrical perspective lambda must be nonzero"))
        iszero(mu_value + lambda_value) && throw(ArgumentError("cylindrical perspective mu + lambda must be nonzero"))
        return new{T}(mu_value, lambda_value)
    end
end

CylindricalPerspective(mu::Real=1, lambda::Real=1) = CylindricalPerspective{promote_type(typeof(float(mu)), typeof(float(lambda)))}(mu, lambda)

"""The MER/Mercator projection."""
struct MercatorProjection <: AbstractSphericalProjection end

"""The SFL/Sanson-Flamsteed projection."""
struct SansonFlamsteedProjection <: AbstractSphericalProjection end

"""The PAR/parabolic projection."""
struct ParabolicProjection <: AbstractSphericalProjection end

"""The MOL/Mollweide projection."""
struct MollweideProjection <: AbstractSphericalProjection end

"""The PCO/polyconic projection."""
struct PolyconicProjection <: AbstractSphericalProjection end

"""The AIT/Hammer-Aitoff projection."""
struct HammerAitoffProjection <: AbstractSphericalProjection end

"""The ZPN/zenithal polynomial projection with coefficients `PV0, PV1, ...`."""
struct ZenithalPolynomial{N,T<:AbstractFloat} <: AbstractZenithalProjection
    coefficients::NTuple{N,T}
    monotonic_limit::T
    radial_limit::T
end

function _zpn_monotonic_limit(coefficients::NTuple{N,T}) where {N,T}
    N <= 2 && return T(pi)
    derivative(value) = evalpoly(value, ntuple(i -> T(i) * coefficients[i + 1], N - 1))
    derivative(zero(T)) > zero(T) || throw(ArgumentError("higher-order zenithal polynomial must have a positive linear coefficient"))
    lower = zero(T)
    for degree in 1:179
        upper = deg2rad(T(degree))
        derivative(upper) <= zero(T) || (lower = upper; continue)
        for _ in 1:64
            middle = (lower + upper) / 2
            derivative(middle) > zero(T) ? (lower = middle) : (upper = middle)
        end
        return (lower + upper) / 2
    end
    return T(pi)
end

function ZenithalPolynomial(coefficients=(0.0, 1.0))
    values = Tuple(float.(coefficients))
    isempty(values) && throw(ArgumentError("zenithal polynomial requires at least one coefficient"))
    all(isfinite, values) || throw(ArgumentError("zenithal polynomial coefficients must be finite"))
    last_nonzero = findlast(!iszero, values)
    last_nonzero === nothing && throw(ArgumentError("zenithal polynomial must have a nonzero coefficient"))
    values = values[1:last_nonzero]
    T = promote_type(map(typeof, values)...)
    typed = T.(values)
    limit = _zpn_monotonic_limit(typed)
    return ZenithalPolynomial{length(typed),T}(typed, limit, _projection_polynomial(typed, limit))
end

abstract type AbstractConicProjection <: AbstractSphericalProjection end

function _conic_angles(::Type{T}, sigma, delta) where {T}
    sigma_value, delta_value = T(sigma), T(delta)
    all(isfinite, (sigma_value, delta_value)) || throw(ArgumentError("conic projection parameters must be finite"))
    abs(sigma_value) <= 90 && abs(delta_value) <= 90 || throw(ArgumentError("conic projection parameters must be between -90 and 90 degrees"))
    return deg2rad(sigma_value), deg2rad(delta_value)
end

function _conic_perspective_constants(sigma::T, delta::T) where {T}
    cone = sin(sigma)
    abs(cone) > _projection_tolerance(T) || throw(ArgumentError("COP is degenerate for sigma = 0"))
    cos_delta = cos(delta)
    abs(cos_delta) > _projection_tolerance(T) || throw(ArgumentError("COP is degenerate for delta = +/-90 degrees"))
    cot_sigma = cos(sigma) / sin(sigma)
    return sigma, cone, rad2deg(cos_delta * cot_sigma), cos_delta, cot_sigma
end

function _conic_equidistant_constants(sigma::T, delta::T) where {T}
    cone = abs(delta) <= _projection_tolerance(T) ? sin(sigma) : sin(sigma) * sin(delta) / delta
    abs(cone) > _projection_tolerance(T) || throw(ArgumentError("COD conic constant is degenerate"))
    return sigma, cone, rad2deg(cos(delta) * cos(sigma) / cone)
end

function _conic_equal_area_constants(sigma::T, delta::T) where {T}
    first_parallel, second_parallel = sigma - delta, sigma + delta
    cone = (sin(first_parallel) + sin(second_parallel)) / 2
    abs(cone) > _projection_tolerance(T) || throw(ArgumentError("COE conic constant is degenerate"))
    scale = rad2deg(one(T)) / cone
    shape = one(T) + sin(first_parallel) * sin(second_parallel)
    y0 = scale * sqrt(max(zero(T), shape - 2cone * sin(sigma)))
    return sigma, cone, scale, shape, y0
end

function _conic_orthomorphic_constants(sigma::T, delta::T) where {T}
    first_parallel, second_parallel = sigma - delta, sigma + delta
    cone = if abs(delta) <= _projection_tolerance(T)
        sin(first_parallel)
    else
        first_tangent = tan((T(pi / 2) - first_parallel) / 2)
        second_tangent = tan((T(pi / 2) - second_parallel) / 2)
        log(cos(second_parallel) / cos(first_parallel)) / log(second_tangent / first_tangent)
    end
    abs(cone) > _projection_tolerance(T) || throw(ArgumentError("COO conic constant is degenerate"))
    first_tangent = tan((T(pi / 2) - first_parallel) / 2)
    scale = rad2deg(cos(first_parallel) / (cone * first_tangent^cone))
    y0 = scale * tan((T(pi / 2) - sigma) / 2)^cone
    return sigma, cone, scale, y0
end

# `constants` holds values precomputed from sigma and delta (in degrees) so evaluation avoids recomputing them.
struct ConicPerspective{T<:AbstractFloat} <: AbstractConicProjection
    sigma::T
    delta::T
    constants::NTuple{5,T}
    ConicPerspective{T}(sigma, delta) where {T} = new{T}(sigma, delta, _conic_perspective_constants(_conic_angles(T, sigma, delta)...))
end

struct ConicEquidistant{T<:AbstractFloat} <: AbstractConicProjection
    sigma::T
    delta::T
    constants::NTuple{3,T}
    ConicEquidistant{T}(sigma, delta) where {T} = new{T}(sigma, delta, _conic_equidistant_constants(_conic_angles(T, sigma, delta)...))
end

struct ConicEqualArea{T<:AbstractFloat} <: AbstractConicProjection
    sigma::T
    delta::T
    constants::NTuple{5,T}
    ConicEqualArea{T}(sigma, delta) where {T} = new{T}(sigma, delta, _conic_equal_area_constants(_conic_angles(T, sigma, delta)...))
end

struct ConicOrthomorphic{T<:AbstractFloat} <: AbstractConicProjection
    sigma::T
    delta::T
    constants::NTuple{4,T}
    ConicOrthomorphic{T}(sigma, delta) where {T} = new{T}(sigma, delta, _conic_orthomorphic_constants(_conic_angles(T, sigma, delta)...))
end

_conic_type(sigma, delta) = promote_type(typeof(float(sigma)), typeof(float(delta)))

ConicPerspective(sigma::Real, delta::Real=0) = ConicPerspective{_conic_type(sigma, delta)}(sigma, delta)
ConicEquidistant(sigma::Real, delta::Real=0) = ConicEquidistant{_conic_type(sigma, delta)}(sigma, delta)
ConicEqualArea(sigma::Real, delta::Real=0) = ConicEqualArea{_conic_type(sigma, delta)}(sigma, delta)
ConicOrthomorphic(sigma::Real, delta::Real=0) = ConicOrthomorphic{_conic_type(sigma, delta)}(sigma, delta)

"""The BON/Bonne projection with standard parallel `theta1`."""
struct BonneProjection{T<:AbstractFloat} <: AbstractSphericalProjection
    theta1::T
    y0::T
    function BonneProjection{T}(theta1) where {T}
        value = T(theta1)
        isfinite(value) && abs(value) <= 90 || throw(ArgumentError("Bonne standard parallel must be between -90 and 90 degrees"))
        # theta1 = 0 is evaluated as SFL, so y0 is unused there.
        radians = deg2rad(value)
        y0 = iszero(value) ? zero(T) : rad2deg(cos(radians) / sin(radians) + radians)
        return new{T}(value, y0)
    end
end

BonneProjection(theta1::Real) = BonneProjection{typeof(float(theta1))}(theta1)

struct TangentialSphericalCube <: AbstractSphericalProjection end
struct COBESphericalCube <: AbstractSphericalProjection end
struct QuadrilateralizedSphericalCube <: AbstractSphericalProjection end

struct HEALPixProjection <: AbstractSphericalProjection
    H::Int
    K::Int
    function HEALPixProjection(H=4, K=3)
        H isa Integer && K isa Integer && H > 0 && K > 0 || throw(ArgumentError("HEALPix H and K must be positive integers"))
        return new(Int(H), Int(K))
    end
end

struct HEALPixPolarProjection <: AbstractSphericalProjection end

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
    # Non-finite coordinates propagate as NaN instead of tripping projection domain checks.
    all(isfinite, coords) || return SVector(T(NaN), T(NaN))
    phi, theta = _intermediate_to_native(transform.projection, T(coords[1]), T(coords[2]))
    # An undefined point is NaN in every component, even when only one formula failed.
    (isnan(phi) || isnan(theta)) && return SVector(oftype(phi, NaN), oftype(theta, NaN))
    return _native_bounds_check(SVector(rad2deg(phi), rad2deg(theta)), coords)
end

# Like WCSLIB's prjbchk, reject native coordinates outside phi in [-180, 180] and theta in
# [-90, 90]; otherwise points beyond a projection's boundary map to wrapped or invalid sky positions.
@inline function _native_bounds_check(native, coords)
    phi, theta = native
    tolerance = max(oftype(phi, 1e-11), 8 * eps(oftype(phi, 180)))
    abs(phi) <= 180 + tolerance && abs(theta) <= 90 + tolerance || return SVector(oftype(phi, NaN), oftype(theta, NaN))
    return SVector(clamp(phi, -oftype(phi, 180), oftype(phi, 180)), clamp(theta, -oftype(theta, 90), oftype(theta, 90)))
end

@inline function (transform::NativeToIntermediate)(coords::StaticVector{2})
    T = _evaluation_type(eltype(coords))
    all(isfinite, coords) || return SVector(T(NaN), T(NaN))
    phi = deg2rad(T(coords[1]))
    theta = deg2rad(T(coords[2]))
    x, y = _native_to_intermediate(transform.projection, phi, theta)
    isnan(x) || isnan(y) ? SVector(oftype(x, NaN), oftype(y, NaN)) : SVector(x, y)
end

# pi/2 - theta with pi/2 in two parts, so a small colatitude keeps its relative accuracy and agrees with cos(theta),
# which the sky rotation uses; the differentiated projection and rotation then agree near the pole.
_half_pi_low(::Type{Float64}) = 6.123233995736766e-17
_half_pi_low(::Type{Float32}) = -4.371139f-8
_half_pi_low(::Type{<:ForwardDiff.Dual{<:Any,V}}) where {V} = _half_pi_low(V)
_half_pi_low(::Type{T}) where {T} = zero(T)
@inline _colatitude(theta) = (oftype(theta, pi / 2) - theta) + _half_pi_low(typeof(theta))

@inline function _zenithal_phi(x, y)
    iszero(x) && iszero(y) && return zero(promote_type(typeof(x), typeof(y)))
    return atan(x, -y)
end

@inline function _zenithal_xy(radius, phi)
    sine, cosine = sincos(phi)
    return radius * sine, -radius * cosine
end

# Coordinates with no image under a projection evaluate to NaN; WCS-level evaluation decides whether to throw or fill.
@inline function _nan_pair(a, b)
    T = promote_type(typeof(a), typeof(b))
    return T(NaN), T(NaN)
end

@inline function _projection_tolerance(::Type{T}) where {T}
    return T(100) * eps(T)
end

function _native_latitude(theta, projection)
    limit = oftype(theta, pi / 2)
    tolerance = _projection_tolerance(typeof(theta))
    abs(theta) <= limit + tolerance || return oftype(theta, NaN)
    return clamp(theta, -limit, limit)
end

@inline _wrap_native_longitude(phi) = phi - oftype(phi, 2pi) * round(phi / oftype(phi, 2pi))

@inline function _intermediate_to_native(::PlateCarreeProjection, x, y)
    abs(y) <= oftype(y, 90) + _projection_tolerance(typeof(y)) || return _nan_pair(x, y)
    return deg2rad(x), deg2rad(clamp(y, oftype(y, -90), oftype(y, 90)))
end

@inline function _native_to_intermediate(projection::PlateCarreeProjection, phi, theta)
    theta = _native_latitude(theta, projection)
    return rad2deg(_wrap_native_longitude(phi)), rad2deg(theta)
end

function _intermediate_to_native(projection::CylindricalEqualArea, x, y)
    T = promote_type(typeof(x), typeof(y))
    argument = T(projection.lambda) * deg2rad(y)
    abs(argument) <= one(T) + _projection_tolerance(T) || return _nan_pair(x, y)
    return deg2rad(x), asin(clamp(argument, -one(T), one(T)))
end

function _native_to_intermediate(projection::CylindricalEqualArea, phi, theta)
    theta = _native_latitude(theta, projection)
    return rad2deg(_wrap_native_longitude(phi)), rad2deg(sin(theta) / typeof(theta)(projection.lambda))
end

function _intermediate_to_native(projection::CylindricalPerspective, x, y)
    T = promote_type(typeof(x), typeof(y))
    mu, lambda = T(projection.mu), T(projection.lambda)
    eta = deg2rad(y) / (mu + lambda)
    argument = eta * mu / hypot(eta, one(T))
    abs(argument) <= one(T) + _projection_tolerance(T) || return _nan_pair(x, y)
    return deg2rad(x) / lambda, atan(eta) + asin(clamp(argument, -one(T), one(T)))
end

function _native_to_intermediate(projection::CylindricalPerspective, phi, theta)
    theta = _native_latitude(theta, projection)
    T = promote_type(typeof(phi), typeof(theta))
    mu, lambda = T(projection.mu), T(projection.lambda)
    denominator = mu + cos(theta)
    iszero(denominator) && return _nan_pair(phi, theta)
    return rad2deg(lambda * _wrap_native_longitude(phi)), rad2deg((mu + lambda) * sin(theta) / denominator)
end

@inline function _intermediate_to_native(::MercatorProjection, x, y)
    T = promote_type(typeof(x), typeof(y))
    return deg2rad(x), 2atan(exp(deg2rad(y))) - oftype(one(T), pi / 2)
end

function _native_to_intermediate(projection::MercatorProjection, phi, theta)
    theta = _native_latitude(theta, projection)
    abs(theta) < oftype(theta, pi / 2) - _projection_tolerance(typeof(theta)) || return _nan_pair(phi, theta)
    return rad2deg(_wrap_native_longitude(phi)), rad2deg(log(tan(oftype(theta, pi / 4) + theta / 2)))
end

function _intermediate_to_native(::SansonFlamsteedProjection, x, y)
    T = promote_type(typeof(x), typeof(y))
    abs(y) <= T(90) + _projection_tolerance(T) || return _nan_pair(x, y)
    theta = deg2rad(clamp(y, T(-90), T(90)))
    cosine = cos(theta)
    if abs(cosine) <= _projection_tolerance(T)
        abs(x) <= _projection_tolerance(T) && return zero(T), theta
        return _nan_pair(x, y)
    end
    return deg2rad(x) / cosine, theta
end

function _native_to_intermediate(projection::SansonFlamsteedProjection, phi, theta)
    theta = _native_latitude(theta, projection)
    return rad2deg(_wrap_native_longitude(phi) * cos(theta)), rad2deg(theta)
end

function _intermediate_to_native(::ParabolicProjection, x, y)
    T = promote_type(typeof(x), typeof(y))
    argument = deg2rad(y) / T(pi)
    abs(argument) <= one(T) + _projection_tolerance(T) || return _nan_pair(x, y)
    theta = 3asin(clamp(argument, -one(T), one(T)))
    abs(theta) <= T(pi / 2) + _projection_tolerance(T) || return _nan_pair(x, y)
    scale = 2cos(2theta / 3) - one(T)
    if abs(scale) <= _projection_tolerance(T)
        abs(x) <= _projection_tolerance(T) && return zero(T), clamp(theta, T(-pi / 2), T(pi / 2))
        return _nan_pair(x, y)
    end
    return deg2rad(x) / scale, theta
end

function _native_to_intermediate(projection::ParabolicProjection, phi, theta)
    theta = _native_latitude(theta, projection)
    return rad2deg(_wrap_native_longitude(phi) * (2cos(2theta / 3) - one(theta))), rad2deg(oftype(theta, pi) * sin(theta / 3))
end

function _mollweide_angle(theta)
    T = typeof(theta)
    isnan(theta) && return theta
    abs(abs(theta) - T(pi / 2)) <= _projection_tolerance(T) && return copysign(T(pi / 2), theta)
    angle = theta
    target = T(pi) * sin(theta)
    for _ in 1:64
        value = 2angle + sin(2angle) - target
        derivative = 2 + 2cos(2angle)
        iszero(derivative) && break
        step = value / derivative
        angle -= step
        abs(step) <= _projection_tolerance(T) && return angle
    end
    error("MOL auxiliary-angle solve failed to converge at native latitude $(rad2deg(theta)) degrees")
end

function _intermediate_to_native(::MollweideProjection, x, y)
    T = promote_type(typeof(x), typeof(y))
    root_two = sqrt(T(2))
    sine = deg2rad(y) / root_two
    abs(sine) <= one(T) + _projection_tolerance(T) || return _nan_pair(x, y)
    angle = asin(clamp(sine, -one(T), one(T)))
    theta = asin(clamp((2angle + sin(2angle)) / T(pi), -one(T), one(T)))
    cosine = cos(angle)
    if abs(cosine) <= _projection_tolerance(T)
        abs(x) <= _projection_tolerance(T) && return zero(T), theta
        return _nan_pair(x, y)
    end
    return deg2rad(x) * T(pi) / (2root_two * cosine), theta
end

function _native_to_intermediate(projection::MollweideProjection, phi, theta)
    theta = _native_latitude(theta, projection)
    T = promote_type(typeof(phi), typeof(theta))
    root_two = sqrt(T(2))
    angle = _mollweide_angle(theta)
    return rad2deg(2root_two / T(pi) * _wrap_native_longitude(phi) * cos(angle)), rad2deg(root_two * sin(angle))
end

function _intermediate_to_native(::PolyconicProjection, x, y)
    T = promote_type(typeof(x), typeof(y))
    tolerance = _projection_tolerance(T)
    abs(y) <= tolerance && return deg2rad(x), zero(T)
    abs(abs(y) - T(90)) <= tolerance && return zero(T), copysign(T(pi / 2), y)
    # |y| > 90 is valid away from the central meridian; the native bounds check rejects invalid points.
    abs(x) <= tolerance && return zero(T), deg2rad(y)
    positive, negative = T(y), zero(T)
    x_squared = T(x)^2
    f_positive, f_negative = x_squared, -x_squared
    theta_degrees = T(y)
    difference = zero(T)
    tangent = zero(T)
    for _ in 1:64
        weight = clamp(f_positive / (f_positive - f_negative), T(0.1), T(0.9))
        theta_degrees = positive - weight * (positive - negative)
        difference = T(y) - theta_degrees
        tangent = tan(deg2rad(theta_degrees))
        value = x_squared + difference * (difference - 2rad2deg(one(T)) / tangent)
        (abs(value) < tolerance || abs(positive - negative) < tolerance) && break
        if value > 0
            positive, f_positive = theta_degrees, value
        else
            negative, f_negative = theta_degrees, value
        end
    end
    first = rad2deg(one(T)) - difference * tangent
    second = T(x) * tangent
    phi_degrees = rad2deg(atan(second, first)) / sin(deg2rad(theta_degrees))
    return deg2rad(phi_degrees), deg2rad(theta_degrees)
end

function _native_to_intermediate(projection::PolyconicProjection, phi, theta)
    theta = _native_latitude(theta, projection)
    wrapped = _wrap_native_longitude(phi)
    abs(theta) <= _projection_tolerance(typeof(theta)) && return rad2deg(wrapped), zero(theta)
    sine = sin(theta)
    cotangent = cos(theta) / sine
    angle = wrapped * sine
    return rad2deg(cotangent * sin(angle)), rad2deg(theta + cotangent * (1 - cos(angle)))
end

@inline function _cos_minus_unnormalized_sinc(value)
    if abs(value) < oftype(value, 0.1)
        # The direct expression loses the leading quadratic term near zero.
        squared = value^2
        series = one(value) / 45360 - squared / 3991680
        series = -one(value) / 840 + squared * series
        series = one(value) / 30 + squared * series
        return squared * (-one(value) / 3 + squared * series)
    end
    return cos(value) - sin(value) / value
end

function transform_jacobian(transform::NativeToIntermediate{<:PolyconicProjection}, coords::StaticVector{2})
    T = _evaluation_type(eltype(coords))
    phi = _wrap_native_longitude(deg2rad(T(coords[1])))
    theta = _native_latitude(deg2rad(T(coords[2])), transform.projection)
    if iszero(theta)
        return SMatrix{2,2}(one(T), zero(T), zero(T), one(T) + phi^2 / 2)
    end

    sine_theta, cosine_theta = sincosd(rad2deg(theta))
    angle = phi * sine_theta
    sine_angle, cosine_angle = sincos(angle)
    inverse_sine = inv(sine_theta)
    longitude_latitude = phi * inverse_sine * (_cos_minus_unnormalized_sinc(angle) - sine_theta^2 * cosine_angle)
    latitude_latitude = one(T) - 2sin(angle / 2)^2 * inverse_sine^2 + phi * cosine_theta^2 * sine_angle * inverse_sine
    return SMatrix{2,2}(cosine_theta * cosine_angle, cosine_theta * sine_angle, longitude_latitude, latitude_latitude)
end

function transform_jacobian(transform::IntermediateToNative{<:PolyconicProjection}, coords::StaticVector{2})
    native = transform(coords)
    return inv(transform_jacobian(NativeToIntermediate(transform.projection), native))
end

function _intermediate_to_native(::HammerAitoffProjection, x, y)
    T = promote_type(typeof(x), typeof(y))
    u, v = deg2rad(T(x)), deg2rad(T(y))
    discriminant = one(T) - u^2 / 8 - v^2 / 2
    discriminant >= -_projection_tolerance(T) || return _nan_pair(x, y)
    z = sqrt(max(discriminant, zero(T)))
    scale = sqrt(2 / (1 + z^2))
    theta = asin(clamp(v / scale, -one(T), one(T)))
    cosine = cos(theta)
    abs(cosine) <= _projection_tolerance(T) && return zero(T), theta
    phi = 2atan(u / (2scale * cosine), z^2 / cosine)
    return phi, theta
end

function _native_to_intermediate(projection::HammerAitoffProjection, phi, theta)
    theta = _native_latitude(theta, projection)
    wrapped = _wrap_native_longitude(phi)
    denominator = one(theta) + cos(theta) * cos(wrapped / 2)
    denominator > zero(theta) || return _nan_pair(phi, theta)
    scale = sqrt(2 / denominator)
    return rad2deg(2scale * cos(theta) * sin(wrapped / 2)), rad2deg(scale * sin(theta))
end

@inline function _projection_polynomial(coefficients, value)
    result = zero(value)
    for coefficient in reverse(coefficients)
        result = muladd(result, value, typeof(value)(coefficient))
    end
    return result
end


function _intermediate_to_native(projection::ZenithalPolynomial{N}, x, y) where {N}
    T = promote_type(typeof(x), typeof(y))
    target = deg2rad(hypot(x, y))
    coefficients = projection.coefficients
    polynomial(value) = _projection_polynomial(coefficients, value) - target
    if N == 1
        abs(polynomial(zero(T))) <= _projection_tolerance(T) || return _nan_pair(x, y)
        distance = zero(T)
    elseif N == 2
        slope = T(coefficients[2])
        distance = (target - T(coefficients[1])) / slope
    else
        lower = zero(T)
        lower_value = polynomial(lower)
        abs(lower_value) <= _projection_tolerance(T) && return _zenithal_phi(x, y), T(pi / 2)
        upper = T(projection.monotonic_limit)
        upper_value = polynomial(upper)
        target < T(coefficients[1]) - _projection_tolerance(T) && return _nan_pair(x, y)
        target > T(projection.radial_limit) + _projection_tolerance(T) && return _nan_pair(x, y)
        for _ in 1:64
            middle = (lower + upper) / 2
            middle_value = polynomial(middle)
            if abs(middle_value) <= _projection_tolerance(T)
                lower = upper = middle
                break
            elseif signbit(middle_value) == signbit(lower_value)
                lower, lower_value = middle, middle_value
            else
                upper, upper_value = middle, middle_value
            end
        end
        distance = (lower + upper) / 2
    end
    zero(T) <= distance <= T(pi) || return _nan_pair(x, y)
    return _zenithal_phi(x, y), T(pi / 2) - distance
end

function _native_to_intermediate(projection::ZenithalPolynomial, phi, theta)
    theta = _native_latitude(theta, projection)
    distance = _colatitude(theta)
    distance <= projection.monotonic_limit + _projection_tolerance(typeof(distance)) || return _nan_pair(phi, theta)
    radius = rad2deg(_projection_polynomial(projection.coefficients, distance))
    return _zenithal_xy(radius, phi)
end

# The native longitude is undefined at a zenithal projection's reference point, the native pole, so the chain rule
# through native coordinates gives NaN there although the projection and sky rotation compose smoothly, and near it
# the projection's Jacobian, taken at the plane point, disagrees with the rotation's, taken at the rounded native
# latitude. Differentiating the projection's native-to-plane direction and the rotation at the same native point keeps
# the two consistent: J = J_rotation * inv(J_plane). At the pole both longitude columns vanish, and derivatives along
# the meridians at native longitudes 0 and 90 take their place.
_is_native_pair(::Type, ::Type) = false
_is_native_pair(::Type{<:IntermediateToNative{<:AbstractZenithalProjection}}, ::Type{<:NativeToCelestial}) = true
_is_native_pair(::Type{<:CelestialToNative}, ::Type{<:NativeToIntermediate{<:AbstractZenithalProjection}}) = true

# Away from the pole along a meridian, so the derivative stays inside the latitude range; the sign cancels in J.
@inline _meridian_derivative(transform, phi::T) where {T} = ForwardDiff.derivative(colatitude -> transform(SVector(phi, T(90) - colatitude)), zero(T))

@inline function _native_derivatives(plane::NativeToIntermediate, sky::NativeToCelestial, native::StaticVector{2,T}) where {T}
    native[2] == 90 || return _transform_jacobian_checked(plane, native), _transform_jacobian_checked(sky, native)
    meridians = (zero(T), T(90))
    return hcat(map(phi -> _meridian_derivative(plane, phi), meridians)...), hcat(map(phi -> _meridian_derivative(sky, phi), meridians)...)
end

@inline function _native_pair_jacobian(projection::IntermediateToNative, rotation::NativeToCelestial, native)
    plane, sky = _native_derivatives(NativeToIntermediate(projection.projection), rotation, native)
    return sky / plane
end

@inline function _native_pair_jacobian(rotation::CelestialToNative, projection::NativeToIntermediate, native)
    plane, sky = _native_derivatives(projection, inverse_transform(rotation), native)
    return plane / sky
end

function transform_jacobian(transform::IntermediateToNative{<:ZenithalPolynomial}, coords::StaticVector{2})
    native = transform(coords)
    inverse_jacobian = _forwarddiff_jacobian(NativeToIntermediate(transform.projection), native)
    return inv(inverse_jacobian)
end

@inline function _conic_radius_angle(x, y, y0, sigma)
    delta_y = y0 - y
    radius = copysign(hypot(x, delta_y), sigma)
    # WCSLIB evaluates atan2(x/r, dy/r) with the signed radius, flipping the angle for southern cones.
    return radius, signbit(sigma) ? atan(-x, -delta_y) : atan(x, delta_y)
end

function _intermediate_to_native(projection::ConicPerspective, x, y)
    T = promote_type(typeof(x), typeof(y))
    sigma, cone, y0, cos_delta, cot_sigma = map(T, projection.constants)
    radius, angle = _conic_radius_angle(T(x), T(y), y0, sigma)
    return angle / cone, sigma + atan(cot_sigma - deg2rad(radius) / cos_delta)
end

function _native_to_intermediate(projection::ConicPerspective, phi, theta)
    theta = _native_latitude(theta, projection)
    T = promote_type(typeof(phi), typeof(theta))
    sigma, cone, y0, cos_delta, _ = map(T, projection.constants)
    angle = cone * _wrap_native_longitude(phi)
    radius = y0 - rad2deg(cos_delta * tan(theta - sigma))
    return radius * sin(angle), -radius * cos(angle) + y0
end

function _intermediate_to_native(projection::ConicEquidistant, x, y)
    T = promote_type(typeof(x), typeof(y))
    sigma, cone, y0 = map(T, projection.constants)
    radius, angle = _conic_radius_angle(T(x), T(y), y0, sigma)
    return angle / cone, deg2rad(y0 + rad2deg(sigma) - radius)
end

function _native_to_intermediate(projection::ConicEquidistant, phi, theta)
    theta = _native_latitude(theta, projection)
    T = promote_type(typeof(phi), typeof(theta))
    sigma, cone, y0 = map(T, projection.constants)
    angle = cone * _wrap_native_longitude(phi)
    radius = y0 + rad2deg(sigma) - rad2deg(theta)
    return radius * sin(angle), -radius * cos(angle) + y0
end

function _intermediate_to_native(projection::ConicEqualArea, x, y)
    T = promote_type(typeof(x), typeof(y))
    sigma, cone, scale, shape, y0 = map(T, projection.constants)
    radius, angle = _conic_radius_angle(T(x), T(y), y0, sigma)
    argument = (scale^2 * shape - radius^2) * cone / (2rad2deg(one(T))^2)
    abs(argument) <= one(T) + _projection_tolerance(T) || return _nan_pair(x, y)
    return angle / cone, asin(clamp(argument, -one(T), one(T)))
end

function _native_to_intermediate(projection::ConicEqualArea, phi, theta)
    theta = _native_latitude(theta, projection)
    T = promote_type(typeof(phi), typeof(theta))
    _, cone, scale, shape, y0 = map(T, projection.constants)
    radius = scale * sqrt(max(zero(T), shape - 2cone * sin(theta)))
    angle = cone * _wrap_native_longitude(phi)
    return radius * sin(angle), -radius * cos(angle) + y0
end

function _intermediate_to_native(projection::ConicOrthomorphic, x, y)
    T = promote_type(typeof(x), typeof(y))
    sigma, cone, scale, y0 = map(T, projection.constants)
    radius, angle = _conic_radius_angle(T(x), T(y), y0, sigma)
    abs(radius) > _projection_tolerance(T) || return _nan_pair(x, y)
    return angle / cone, T(pi / 2) - 2atan((radius / scale)^(inv(cone)))
end

function _native_to_intermediate(projection::ConicOrthomorphic, phi, theta)
    theta = _native_latitude(theta, projection)
    T = promote_type(typeof(phi), typeof(theta))
    _, cone, scale, y0 = map(T, projection.constants)
    theta > T(-pi / 2) || return _nan_pair(phi, theta)
    radius = scale * tan((T(pi / 2) - theta) / 2)^cone
    angle = cone * _wrap_native_longitude(phi)
    return radius * sin(angle), -radius * cos(angle) + y0
end

function _intermediate_to_native(projection::BonneProjection, x, y)
    T = promote_type(typeof(x), typeof(y))
    theta1 = T(projection.theta1)
    iszero(theta1) && return _intermediate_to_native(SansonFlamsteedProjection(), x, y)
    y0 = T(projection.y0)
    radius, angle = _conic_radius_angle(T(x), T(y), y0, theta1)
    theta = deg2rad(y0 - radius)
    abs(theta) <= T(pi / 2) + _projection_tolerance(T) || return _nan_pair(x, y)
    cosine = cos(theta)
    phi = abs(cosine) <= _projection_tolerance(T) ? zero(T) : angle * deg2rad(radius) / cosine
    return phi, clamp(theta, T(-pi / 2), T(pi / 2))
end

function _native_to_intermediate(projection::BonneProjection, phi, theta)
    theta = _native_latitude(theta, projection)
    T = promote_type(typeof(phi), typeof(theta))
    theta1 = T(projection.theta1)
    iszero(theta1) && return _native_to_intermediate(SansonFlamsteedProjection(), phi, theta)
    y0 = T(projection.y0)
    radius = y0 - rad2deg(theta)
    angle = abs(radius) <= _projection_tolerance(T) ? zero(T) : _wrap_native_longitude(phi) * cos(theta) * rad2deg(one(T)) / radius
    return radius * sin(angle), -radius * cos(angle) + y0
end

const _CUBE_FACE_X = (0, 0, 2, 4, 6, 0)
const _CUBE_FACE_Y = (2, 0, 0, 0, 0, -2)

@inline function _sphere_vector(phi, theta)
    cosine = cos(theta)
    return cosine * cos(phi), cosine * sin(phi), sin(theta)
end

function _cube_face(longitude, transverse, vertical)
    face, dominant = 0, vertical
    longitude > dominant && ((face, dominant) = (1, longitude))
    transverse > dominant && ((face, dominant) = (2, transverse))
    -longitude > dominant && ((face, dominant) = (3, -longitude))
    -transverse > dominant && ((face, dominant) = (4, -transverse))
    -vertical > dominant && ((face, dominant) = (5, -vertical))
    if face == 1
        x, y = transverse / dominant, vertical / dominant
    elseif face == 2
        x, y = -longitude / dominant, vertical / dominant
    elseif face == 3
        x, y = -transverse / dominant, vertical / dominant
    elseif face == 4
        x, y = longitude / dominant, vertical / dominant
    elseif face == 5
        x, y = transverse / dominant, longitude / dominant
    else
        x, y = transverse / dominant, -longitude / dominant
    end
    return face, x, y
end

function _cube_canvas_face(x, y)
    T = promote_type(typeof(x), typeof(y))
    face_x, face_y = T(x) / 45, T(y) / 45
    face_x < -one(T) && (face_x += T(8))
    tolerance = _projection_tolerance(T)
    if face_y > one(T) + tolerance
        face, face_y = 0, face_y - 2
    elseif face_y < -one(T) - tolerance
        face, face_y = 5, face_y + 2
    elseif face_x > 5 + tolerance
        face, face_x = 4, face_x - 6
    elseif face_x > 3 + tolerance
        face, face_x = 3, face_x - 4
    elseif face_x > 1 + tolerance
        face, face_x = 2, face_x - 2
    else
        face = 1
    end
    max(abs(face_x), abs(face_y)) <= one(T) + tolerance || return 1, oftype(face_x, NaN), oftype(face_y, NaN)
    return face, clamp(face_x, -one(T), one(T)), clamp(face_y, -one(T), one(T))
end

function _cube_vector(face, x, y)
    if face == 1
        longitude, transverse, vertical = one(x), x, y
    elseif face == 2
        longitude, transverse, vertical = -x, one(x), y
    elseif face == 3
        longitude, transverse, vertical = -one(x), -x, y
    elseif face == 4
        longitude, transverse, vertical = x, -one(x), y
    elseif face == 5
        longitude, transverse, vertical = y, x, -one(x)
    else
        longitude, transverse, vertical = -y, x, one(x)
    end
    norm = hypot(longitude, transverse, vertical)
    return longitude / norm, transverse / norm, vertical / norm
end

@inline function _cube_native(vector)
    longitude, transverse, vertical = vector
    return atan(transverse, longitude), asin(clamp(vertical, -one(vertical), one(vertical)))
end

function _native_to_intermediate(projection::TangentialSphericalCube, phi, theta)
    theta = _native_latitude(theta, projection)
    face, x, y = _cube_face(_sphere_vector(phi, theta)...)
    return 45 * (x + _CUBE_FACE_X[face + 1]), 45 * (y + _CUBE_FACE_Y[face + 1])
end

function _intermediate_to_native(::TangentialSphericalCube, x, y)
    face, face_x, face_y = _cube_canvas_face(x, y)
    return _cube_native(_cube_vector(face, face_x, face_y))
end

function _cobe_forward(primary, secondary)
    T = promote_type(typeof(primary), typeof(secondary))
    primary_squared, secondary_squared = primary^2, secondary^2
    primary_complement, secondary_complement = one(T) - primary_squared, one(T) - secondary_squared
    product = abs(primary * secondary) > T(1e-16) ? primary_squared * secondary_squared : zero(T)
    primary_fourth = primary_squared > T(1e-16) ? primary_squared^2 : zero(T)
    secondary_fourth = secondary_squared > T(1e-16) ? secondary_squared^2 : zero(T)
    return primary * (primary_squared + primary_complement * (T(1.37484847732) + secondary_squared *
        (T(-0.13161671474) * primary_complement + T(0.004869491981) * primary_squared + secondary_complement *
        (T(0.141189631152) + T(0.0809701286525) * primary_squared - T(0.281528535557) * secondary_squared +
        T(0.15384112876) * product - T(0.178251207466) * primary_fourth + T(0.106959469314) * secondary_fourth)) +
        primary_squared * (T(-0.159596235474) - primary_complement * (T(0.0759196200467) - T(0.0217762490699) * primary_squared))))
end

function _cobe_inverse(projected, secondary)
    T = promote_type(typeof(projected), typeof(secondary))
    x, y = projected^2, secondary^2
    z0 = T(-0.27292696) + x * (T(-0.07629969) + x * (T(-0.22797056) + x * (T(0.54852384) + x * (T(-0.62930065) + x * (T(0.25795794) + x * T(0.02584375))))))
    z1 = T(-0.02819452) + x * (T(-0.01471565) + x * (T(0.48051509) + x * (T(-1.74114454) + x * (T(1.71547508) - x * T(0.53022337)))))
    z2 = T(0.27058160) + x * (T(-0.56800938) + x * (T(0.30803317) + x * (T(0.98938102) - x * T(0.83180469))))
    z3 = T(-0.60441560) + x * (T(1.50880086) + x * (T(-0.93678576) + x * T(0.08693841)))
    z4 = T(0.93412077) + x * (T(-1.41601920) + x * T(0.33887446))
    z5 = T(-0.63915306) + x * T(0.52032238)
    correction = z0 + y * (z1 + y * (z2 + y * (z3 + y * (z4 + y * (z5 + y * T(0.14381585))))))
    return projected + projected * (one(T) - x) * correction
end

function _native_to_intermediate(projection::COBESphericalCube, phi, theta)
    theta = _native_latitude(theta, projection)
    face, x, y = _cube_face(_sphere_vector(phi, theta)...)
    return 45 * (_cobe_forward(x, y) + _CUBE_FACE_X[face + 1]), 45 * (_cobe_forward(y, x) + _CUBE_FACE_Y[face + 1])
end

function _intermediate_to_native(::COBESphericalCube, x, y)
    face, face_x, face_y = _cube_canvas_face(x, y)
    return _cube_native(_cube_vector(face, _cobe_inverse(face_x, face_y), _cobe_inverse(face_y, face_x)))
end


function _qsc_face_coordinates(first, second, complement)
    T = promote_type(typeof(first), typeof(second), typeof(complement))
    iszero(first) && iszero(second) && return zero(T), zero(T)
    if abs(first) >= abs(second)
        ratio = second / first
        magnitude = sqrt(complement / (one(T) - inv(sqrt(2 + ratio^2))))
        x = copysign(magnitude, first)
        y = x / 15 * (atand(ratio) - asind(ratio / sqrt(2 * (1 + ratio^2))))
    else
        ratio = first / second
        magnitude = sqrt(complement / (one(T) - inv(sqrt(2 + ratio^2))))
        y = copysign(magnitude, second)
        x = y / 15 * (atand(ratio) - asind(ratio / sqrt(2 * (1 + ratio^2))))
    end
    return x, y
end

function _native_to_intermediate(projection::QuadrilateralizedSphericalCube, phi, theta)
    theta = _native_latitude(theta, projection)
    longitude, transverse, vertical = _sphere_vector(phi, theta)
    face, _, _ = _cube_face(longitude, transverse, vertical)
    dominant = max(abs(longitude), abs(transverse), abs(vertical))
    if face == 1
        first, second = transverse, vertical
    elseif face == 2
        first, second = -longitude, vertical
    elseif face == 3
        first, second = -transverse, vertical
    elseif face == 4
        first, second = longitude, vertical
    elseif face == 5
        first, second = transverse, longitude
    else
        first, second = transverse, -longitude
    end
    x, y = _qsc_face_coordinates(first, second, one(dominant) - dominant)
    return 45 * (x + _CUBE_FACE_X[face + 1]), 45 * (y + _CUBE_FACE_Y[face + 1])
end


function _intermediate_to_native(::QuadrilateralizedSphericalCube, x, y)
    face, face_x, face_y = _cube_canvas_face(x, y)
    T = promote_type(typeof(face_x), typeof(face_y))
    if abs(face_x) <= _projection_tolerance(T) && abs(face_y) <= _projection_tolerance(T)
        dominant, first, second = one(T), zero(T), zero(T)
    else
        direct = abs(face_x) > abs(face_y)
        angle = direct ? 15face_y / face_x : 15face_x / face_y
        ratio = sind(angle) / (cosd(angle) - inv(sqrt(T(2))))
        complement = (direct ? face_x^2 : face_y^2) * (one(T) - inv(sqrt(2 + ratio^2)))
        dominant = one(T) - complement
        magnitude = sqrt(max(zero(T), complement * (2 - complement) / (1 + ratio^2)))
        if direct
            first, second = copysign(magnitude, face_x), copysign(magnitude, face_x) * ratio
        else
            second, first = copysign(magnitude, face_y), copysign(magnitude, face_y) * ratio
        end
    end
    vector = face == 1 ? (dominant, first, second) : face == 2 ? (-first, dominant, second) :
        face == 3 ? (-dominant, -first, second) : face == 4 ? (first, -dominant, second) :
        face == 5 ? (second, first, -dominant) : (-second, first, dominant)
    norm = hypot(vector...)
    return _cube_native(map(value -> value / norm, vector))
end

function _native_to_intermediate(projection::HEALPixProjection, phi, theta)
    theta = _native_latitude(theta, projection)
    T = promote_type(typeof(phi), typeof(theta))
    H, K = T(projection.H), T(projection.K)
    sine = sin(theta)
    abs(sine) <= (K - 1) / K && return rad2deg(_wrap_native_longitude(phi)), 90K / H * sine
    sigma = sqrt(K * (1 - abs(sine)))
    longitude = rad2deg(_wrap_native_longitude(phi))
    step = 180 / H
    center = -180 + (2floor((longitude + 180) * H / 360) + 1) * step
    difference = longitude - center
    if iseven(projection.K) && sine < 0
        facet = floor(Int, longitude / step) + (projection.H % 2)
        difference += isodd(facet) ? -step : step
    end
    x = longitude + difference * (sigma - 1)
    y = copysign(step * ((K + 1) / 2 - sigma), sine)
    x > 180 && (x = 360 - x)
    return x, y
end

function _intermediate_to_native(projection::HEALPixProjection, x, y)
    T = promote_type(typeof(x), typeof(y))
    H, K = T(projection.H), T(projection.K)
    boundary, scale, step = 90 * (K - 1) / H, 90K / H, 180 / H
    if abs(y) <= boundary
        return deg2rad(x), asin(clamp(y / scale, -one(T), one(T)))
    end
    sigma = (K + 1) / 2 - abs(y) / step
    zero(T) <= sigma <= one(T) + _projection_tolerance(T) || return _nan_pair(x, y)
    theta = copysign(asin(clamp(1 - sigma^2 / K, -one(T), one(T))), y)
    center = -180 + (2floor((x + 180) * H / 360) + 1) * step
    difference = x - center
    if iseven(projection.K) && y < 0
        facet = floor(Int, x / step) + (projection.H % 2)
        difference += isodd(facet) ? -step : step
    end
    # As in WCSLIB's hpxx2s, points between the polar facet triangles are not on the sphere.
    abs(difference) < (step + T(1e-12)) * sigma || iszero(difference) || return _nan_pair(x, y)
    phi = iszero(sigma) ? zero(T) : deg2rad(x + difference * (inv(sigma) - 1))
    return phi, theta
end

function _native_to_intermediate(projection::HEALPixPolarProjection, phi, theta)
    theta = _native_latitude(theta, projection)
    T = promote_type(typeof(phi), typeof(theta))
    sine = sin(theta)
    longitude = mod(rad2deg(phi) + 180, T(360)) - 180
    shifted = longitude + 180
    within = mod(shifted, T(90))
    if abs(sine) <= T(2 / 3)
        first, second = within, T(67.5) * sine
    else
        sigma = sqrt(3 * (1 - abs(sine)))
        first, second = 45 + (within - 45) * sigma, copysign(45 * (2 - sigma), sine)
    end
    first -= 45
    second -= 90
    scale = inv(sqrt(T(2)))
    if longitude < -90
        return scale * (-first + second), scale * (-first - second)
    elseif longitude < 0
        return scale * (first + second), scale * (-first + second)
    elseif longitude < 90
        return scale * (first - second), scale * (first + second)
    end
    return scale * (-first - second), scale * (first - second)
end

function _intermediate_to_native(::HEALPixPolarProjection, x, y)
    T = promote_type(typeof(x), typeof(y))
    scaled_x, scaled_y = T(x) / sqrt(T(2)), T(y) / sqrt(T(2))
    if scaled_x <= 0 && scaled_y > 0
        first, second, base = -scaled_x - scaled_y, scaled_x - scaled_y, T(-180)
    elseif scaled_x < 0 && scaled_y <= 0
        first, second, base = scaled_x - scaled_y, scaled_x + scaled_y, T(-90)
    elseif scaled_x >= 0 && scaled_y < 0
        first, second, base = scaled_x + scaled_y, -scaled_x + scaled_y, zero(T)
    else
        first, second, base = -scaled_x + scaled_y, -scaled_x - scaled_y, T(90)
    end
    facet_offset = first
    first += 45
    second += 90
    abs_second = abs(second)
    # Facet checks from WCSLIB's xphx2s; points outside the facets are not on the sphere.
    outside = abs_second <= 45 ? abs(facet_offset) > 45 + T(1e-12) : second < -45 && abs(facet_offset) > second + 90 + T(1e-12)
    outside && return _nan_pair(x, y)
    if abs_second <= 45
        return deg2rad(base + first), asin(clamp(second / T(67.5), -one(T), one(T)))
    elseif abs_second <= 90 + _projection_tolerance(T)
        sigma = (90 - abs_second) / 45
        phi = sigma <= _projection_tolerance(T) ? zero(T) : deg2rad(base + 45 + (first - 45) / sigma)
        theta = copysign(asin(clamp(1 - sigma^2 / 3, -one(T), one(T))), second)
        return phi, theta
    end
    return _nan_pair(x, y)
end

@inline function _intermediate_to_native(::GnomonicProjection, x, y)
    T = promote_type(typeof(x), typeof(y))
    return _zenithal_phi(x, y), atan(rad2deg(one(T)), hypot(x, y))
end

@inline function _native_to_intermediate(projection::GnomonicProjection, phi, theta)
    theta = _native_latitude(theta, projection)
    sine, cosine = sincos(theta)
    sine > zero(sine) || return _nan_pair(phi, theta)
    return _zenithal_xy(rad2deg(cosine / sine), phi)
end

@inline function _intermediate_to_native(::StereographicProjection, x, y)
    T = promote_type(typeof(x), typeof(y))
    # WCSLIB's 90 - 2 atan(R / 2r0); the equivalent asin((1 - r^2) / (1 + r^2)) loses the latitude near the pole.
    return _zenithal_phi(x, y), T(pi / 2) - 2atan(deg2rad(hypot(x, y)) / 2)
end

@inline function _native_to_intermediate(projection::StereographicProjection, phi, theta)
    theta = _native_latitude(theta, projection)
    denominator = one(theta) + sin(theta)
    iszero(denominator) && return _nan_pair(phi, theta)
    return _zenithal_xy(rad2deg(2cos(theta) / denominator), phi)
end

@inline function _intermediate_to_native(::ZenithalEquidistant, x, y)
    T = promote_type(typeof(x), typeof(y))
    theta = oftype(one(T), pi / 2) - deg2rad(hypot(x, y))
    theta >= -oftype(theta, pi / 2) - _projection_tolerance(T) || return _nan_pair(x, y)
    return _zenithal_phi(x, y), max(theta, -oftype(theta, pi / 2))
end

@inline function _native_to_intermediate(projection::ZenithalEquidistant, phi, theta)
    theta = _native_latitude(theta, projection)
    return _zenithal_xy(rad2deg(_colatitude(theta)), phi)
end

@inline function _intermediate_to_native(::ZenithalEqualArea, x, y)
    T = promote_type(typeof(x), typeof(y))
    argument = deg2rad(hypot(x, y)) / 2
    tolerance = _projection_tolerance(T)
    argument <= one(T) + tolerance || return _nan_pair(x, y)
    theta = oftype(one(T), pi / 2) - 2asin(clamp(argument, zero(T), one(T)))
    return _zenithal_phi(x, y), theta
end

@inline function _native_to_intermediate(projection::ZenithalEqualArea, phi, theta)
    theta = _native_latitude(theta, projection)
    radius = rad2deg(2sin(_colatitude(theta) / 2))
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
        radius_squared <= one(T) + _projection_tolerance(T) || return _nan_pair(x, y)
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
        discriminant >= -tolerance || return _nan_pair(x, y)
        root = sqrt(max(discriminant, zero(T)))
        sine1 = (-b + root) / a
        sine2 = (-b - root) / a
        sine = max(sine1, sine2)
        if sine > one(T) + tolerance
            sine = min(sine1, sine2)
        end
        -one(T) - tolerance <= sine <= one(T) + tolerance || return _nan_pair(x, y)
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
        theta >= zero(T) || return _nan_pair(phi, theta)
    else
        minimum_theta = -atan(xi * sine_phi - eta * cosine_phi)
        theta >= minimum_theta - _projection_tolerance(T) || return _nan_pair(phi, theta)
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
    iszero(q) && return _nan_pair(x, y)
    ratio = radius / q
    argument = ratio * mu / sqrt(ratio^2 + one(T))
    tolerance = _projection_tolerance(T)
    abs(argument) <= one(T) + tolerance || return _nan_pair(x, y)
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
    iszero(denominator) && return _nan_pair(phi, theta)

    theta_degrees = rad2deg(theta)
    theta_degrees >= T(projection.divergence_latitude) - _projection_tolerance(T) || return _nan_pair(phi, theta)
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
            theta >= max(first_limit, second_limit) - _projection_tolerance(T) || return _nan_pair(phi, theta)
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
        discriminant >= -tolerance || return _nan_pair(x, y)
        root = sqrt(max(discriminant, zero(T)))
        sine1 = (-b + root) / a
        sine2 = (-b - root) / a
        sine = max(sine1, sine2)
        if sine > one(T) + tolerance
            sine = min(sine1, sine2)
        end
        -one(T) - tolerance <= sine <= one(T) + tolerance || return _nan_pair(x, y)
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
    iszero(denominator) && return _nan_pair(phi, theta)

    theta_degrees = rad2deg(theta)
    theta_degrees >= T(projection.divergence_latitude) - _projection_tolerance(T) || return _nan_pair(phi, theta)
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
                theta >= max(first_limit, second_limit) - _projection_tolerance(T) || return _nan_pair(phi, theta)
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
    xi = _colatitude(theta) / 2
    xi >= zero(T) || return oftype(theta, NaN)
    xi < oftype(theta, pi / 2) || return oftype(theta, NaN)
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
    target <= upper_radius || return _nan_pair(x, y)
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
native_fiducial(::Union{PlateCarreeProjection,CylindricalEqualArea,CylindricalPerspective,MercatorProjection,SansonFlamsteedProjection,ParabolicProjection,MollweideProjection,PolyconicProjection,HammerAitoffProjection}) = (0.0, 0.0)
native_fiducial(projection::AbstractConicProjection) = (0.0, projection.sigma)
native_fiducial(::BonneProjection) = (0.0, 0.0)
native_fiducial(::Union{TangentialSphericalCube,COBESphericalCube,QuadrilateralizedSphericalCube,HEALPixProjection}) = (0.0, 0.0)
native_fiducial(::HEALPixPolarProjection) = (0.0, 90.0)

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
    celestial_transform(projection, reference; lonpole=nothing, latpole=90, native_reference=native_fiducial(projection), offset=false)

Construct an intermediate-plane-to-celestial transform using a two-coordinate
celestial reference point and FITS Paper II native-pole conventions. All
coordinates and pole parameters are in degrees. `native_reference` specifies
the native longitude and latitude of the fiducial point. When `offset` is true,
the projection is translated so that this native fiducial point has zero
intermediate coordinates, as specified by FITS longitude-axis `PVi_0`.
"""
function celestial_transform(projection::AbstractSphericalProjection, reference; lonpole=nothing, latpole=90, native_reference=native_fiducial(projection), offset=false)
    reference isa Tuple || reference isa AbstractVector || throw(ArgumentError("celestial reference coordinate must be a tuple or vector"))
    length(reference) == 2 || throw(DimensionMismatch("celestial reference coordinate has length $(length(reference)), expected 2"))
    all(value -> value isa Real, reference) || throw(ArgumentError("celestial reference coordinates must be real numbers"))
    native_reference isa Tuple || native_reference isa AbstractVector || throw(ArgumentError("native reference coordinate must be a tuple or vector"))
    length(native_reference) == 2 || throw(DimensionMismatch("native reference coordinate has length $(length(native_reference)), expected 2"))
    all(value -> value isa Real, native_reference) || throw(ArgumentError("native reference coordinates must be real numbers"))
    lonpole === nothing || lonpole isa Real || throw(ArgumentError("lonpole must be a real number or nothing"))
    latpole isa Real || throw(ArgumentError("latpole must be a real number"))
    offset isa Bool || throw(ArgumentError("offset must be a boolean"))
    phi0, theta0 = native_reference
    reference_longitude, reference_latitude = reference
    phi_p = lonpole === nothing ? default_lonpole(reference_latitude, theta0, phi0) : lonpole
    alpha_p, delta_p = _native_pole(reference_longitude, reference_latitude, phi0, theta0, phi_p, latpole)
    projection_transform = if offset
        plane_reference = NativeToIntermediate(projection)(SVector(phi0, theta0))
        AffineTransform(SMatrix{2,2}(I), plane_reference) | IntermediateToNative(projection)
    else
        IntermediateToNative(projection)
    end
    return projection_transform | NativeToCelestial(alpha_p, delta_p, phi_p)
end
