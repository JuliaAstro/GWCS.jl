const _FITS_SPEED_OF_LIGHT = 2.99792458e8
const _FITS_PLANCK_CONSTANT = 6.62607015e-34

@inline function _fits_refractive_index_air(inverse_wavelength)
    squared = inverse_wavelength^2
    return oftype(inverse_wavelength, 1.000064328) + oftype(inverse_wavelength, 2.554e8) / (oftype(inverse_wavelength, 0.41e14) - squared) +
        oftype(inverse_wavelength, 294.981e8) / (oftype(inverse_wavelength, 1.46e14) - squared)
end

@inline function _fits_wave_to_air(wavelength)
    index = one(wavelength)
    for _ in 1:4
        index = _fits_refractive_index_air(index / wavelength)
    end
    return wavelength / index
end

@inline function _fits_air_to_wave(wavelength)
    return wavelength * _fits_refractive_index_air(inv(wavelength))
end

@inline function _fits_frequency_to_velocity(frequency, rest_frequency)
    ratio = rest_frequency / frequency
    return oftype(frequency, _FITS_SPEED_OF_LIGHT) * (ratio^2 - 1) / (ratio^2 + 1)
end

@inline function _fits_velocity_to_frequency(velocity, rest_frequency)
    light = oftype(velocity, _FITS_SPEED_OF_LIGHT)
    return rest_frequency * sqrt((light - velocity) / (light + velocity))
end

@inline function _fits_wavelength_to_velocity(wavelength, rest_wavelength)
    squared = wavelength^2
    rest_squared = rest_wavelength^2
    return oftype(wavelength, _FITS_SPEED_OF_LIGHT) * (squared - rest_squared) / (squared + rest_squared)
end

@inline function _fits_velocity_to_wavelength(velocity, rest_wavelength)
    light = oftype(velocity, _FITS_SPEED_OF_LIGHT)
    return rest_wavelength * sqrt((light + velocity) / (light - velocity))
end

@inline function _fits_spectral_to_parent(value, ::Val{S}, rest_frequency, rest_wavelength) where {S}
    S === :FREQ && return value
    S === :AFRQ && return value / oftype(value, 2pi)
    S === :ENER && return value / oftype(value, _FITS_PLANCK_CONSTANT)
    S === :WAVN && return value * oftype(value, _FITS_SPEED_OF_LIGHT)
    S === :VRAD && return rest_frequency * (1 - value / oftype(value, _FITS_SPEED_OF_LIGHT))
    S === :WAVE && return value
    S === :VOPT && return rest_wavelength * (1 + value / oftype(value, _FITS_SPEED_OF_LIGHT))
    S === :ZOPT && return rest_wavelength * (1 + value)
    S === :AWAV && return value
    S === :VELO && return value
    S === :BETA && return value * oftype(value, _FITS_SPEED_OF_LIGHT)
    throw(ArgumentError("unsupported FITS spectral coordinate type $S"))
end

@inline function _fits_parent_to_spectral(value, ::Val{S}, rest_frequency, rest_wavelength) where {S}
    S === :FREQ && return value
    S === :AFRQ && return oftype(value, 2pi) * value
    S === :ENER && return oftype(value, _FITS_PLANCK_CONSTANT) * value
    S === :WAVN && return value / oftype(value, _FITS_SPEED_OF_LIGHT)
    S === :VRAD && return oftype(value, _FITS_SPEED_OF_LIGHT) * (rest_frequency - value) / rest_frequency
    S === :WAVE && return value
    S === :VOPT && return oftype(value, _FITS_SPEED_OF_LIGHT) * (value - rest_wavelength) / rest_wavelength
    S === :ZOPT && return (value - rest_wavelength) / rest_wavelength
    S === :AWAV && return value
    S === :VELO && return value
    S === :BETA && return value / oftype(value, _FITS_SPEED_OF_LIGHT)
    throw(ArgumentError("unsupported FITS spectral coordinate type $S"))
end

@inline function _fits_x_to_parent(value, ::Val{X}, ::Val{P}, rest_frequency, rest_wavelength) where {X,P}
    X === P && return value
    X === :F && P === :W && return oftype(value, _FITS_SPEED_OF_LIGHT) / value
    X === :W && P === :F && return oftype(value, _FITS_SPEED_OF_LIGHT) / value
    X === :F && P === :V && return _fits_frequency_to_velocity(value, rest_frequency)
    X === :V && P === :F && return _fits_velocity_to_frequency(value, rest_frequency)
    X === :W && P === :V && return _fits_wavelength_to_velocity(value, rest_wavelength)
    X === :V && P === :W && return _fits_velocity_to_wavelength(value, rest_wavelength)
    X === :F && P === :A && return _fits_wave_to_air(oftype(value, _FITS_SPEED_OF_LIGHT) / value)
    X === :W && P === :A && return _fits_wave_to_air(value)
    X === :V && P === :A && return _fits_wave_to_air(_fits_velocity_to_wavelength(value, rest_wavelength))
    X === :A && P === :F && return oftype(value, _FITS_SPEED_OF_LIGHT) / _fits_air_to_wave(value)
    X === :A && P === :W && return _fits_air_to_wave(value)
    X === :A && P === :V && return _fits_wavelength_to_velocity(_fits_air_to_wave(value), rest_wavelength)
    throw(ArgumentError("unsupported FITS spectral conversion $X to $P"))
end

@inline _fits_parent_to_x(value, x, parent, rest_frequency, rest_wavelength) = _fits_x_to_parent(value, parent, x, rest_frequency, rest_wavelength)

struct FITSLogTransform{T<:AbstractFloat} <: AbstractTransform{1,1}
    reference_value::T
    output_scale::T
end

function FITSLogTransform(reference_value; output_scale=1)
    reference, scale = promote(float(reference_value), float(output_scale))
    isfinite(reference) && reference > 0 || throw(ArgumentError("FITS logarithmic reference value must be finite and positive"))
    isfinite(scale) && scale > 0 || throw(ArgumentError("FITS logarithmic output scale must be finite and positive"))
    return FITSLogTransform{typeof(reference)}(reference, scale)
end

@inline function (transform::FITSLogTransform)(coords::StaticVector{1})
    offset = only(coords)
    reference = oftype(offset, transform.reference_value)
    return SVector(reference * exp(offset / reference) / oftype(offset, transform.output_scale))
end

struct FITSLogInverse{F<:FITSLogTransform} <: AbstractTransform{1,1}
    forward::F
end

@inline function (transform::FITSLogInverse)(coords::StaticVector{1})
    forward = transform.forward
    world = only(coords) * oftype(only(coords), forward.output_scale)
    world > 0 || return SVector(oftype(world, NaN))
    reference = oftype(world, forward.reference_value)
    return SVector(reference * log(world / reference))
end

inverse_transform(transform::FITSLogTransform) = FITSLogInverse(transform)
inverse_transform(transform::FITSLogInverse) = transform.forward

struct FITSSpectralTransform{X,P,S,A,T<:AbstractFloat} <: AbstractTransform{1,1}
    reference_value::T
    x_reference::T
    rest_frequency::T
    rest_wavelength::T
    output_scale::T
end

function FITSSpectralTransform(::Val{X}, ::Val{P}, ::Val{S}, ::Val{A}, reference_value, rest_frequency, rest_wavelength; output_scale=1) where {X,P,S,A}
    values = promote(float(reference_value), float(rest_frequency), float(rest_wavelength), float(output_scale))
    reference, frequency, wavelength, scale = values
    !isfinite(frequency) && isfinite(wavelength) && (frequency = oftype(wavelength, _FITS_SPEED_OF_LIGHT) / wavelength)
    !isfinite(wavelength) && isfinite(frequency) && (wavelength = oftype(frequency, _FITS_SPEED_OF_LIGHT) / frequency)
    parent = _fits_spectral_to_parent(reference, Val(S), frequency, wavelength)
    x_reference = _fits_parent_to_x(parent, Val(X), Val(P), frequency, wavelength)
    all(isfinite, (reference, x_reference, scale)) || throw(ArgumentError("FITS spectral reference values do not define a finite transform"))
    return FITSSpectralTransform{X,P,S,A,typeof(reference)}(reference, x_reference, frequency, wavelength, scale)
end

@inline function (transform::FITSSpectralTransform{X,P,S,A})(coords::StaticVector{1}) where {X,P,S,A}
    offset = only(coords)
    absolute = offset + oftype(offset, transform.x_reference)
    parent = _fits_x_to_parent(absolute, Val(X), Val(P), oftype(offset, transform.rest_frequency), oftype(offset, transform.rest_wavelength))
    world = _fits_parent_to_spectral(parent, Val(S), oftype(offset, transform.rest_frequency), oftype(offset, transform.rest_wavelength))
    return SVector(world / oftype(world, transform.output_scale))
end

struct FITSSpectralInverse{X,P,S,A,F} <: AbstractTransform{1,1}
    forward::F
end

@inline function (transform::FITSSpectralInverse{X,P,S,A})(coords::StaticVector{1}) where {X,P,S,A}
    forward = transform.forward
    world = only(coords) * oftype(only(coords), forward.output_scale)
    parent = _fits_spectral_to_parent(world, Val(S), oftype(world, forward.rest_frequency), oftype(world, forward.rest_wavelength))
    offset = _fits_parent_to_x(parent, Val(X), Val(P), oftype(world, forward.rest_frequency), oftype(world, forward.rest_wavelength)) - oftype(world, forward.x_reference)
    return SVector(offset)
end

inverse_transform(transform::FITSSpectralTransform{X,P,S,A}) where {X,P,S,A} = FITSSpectralInverse{X,P,S,A,typeof(transform)}(transform)
inverse_transform(transform::FITSSpectralInverse) = transform.forward

struct FITSGrismTransform{X,P,S,A,T<:AbstractFloat} <: AbstractTransform{1,1}
    gamma_reference::T
    beta_reference::T
    beta_offset::T
    constant::T
    denominator::T
    rest_frequency::T
    rest_wavelength::T
    output_scale::T
    linear_scale::T
end

function FITSGrismTransform(::Val{X}, ::Val{P}, ::Val{S}, ::Val{A}, reference_value, rest_frequency, rest_wavelength, parameters; output_scale=1) where {X,P,S,A}
    length(parameters) == 7 || throw(DimensionMismatch("FITS grism transform requires seven PV parameters"))
    values = promote(float(reference_value), float(rest_frequency), float(rest_wavelength), float(output_scale), float.(parameters)...)
    reference, frequency, wavelength, scale = values[1:4]
    parameters = values[5:end]
    !isfinite(frequency) && isfinite(wavelength) && (frequency = oftype(wavelength, _FITS_SPEED_OF_LIGHT) / wavelength)
    !isfinite(wavelength) && isfinite(frequency) && (wavelength = oftype(frequency, _FITS_SPEED_OF_LIGHT) / frequency)
    parent_reference = _fits_spectral_to_parent(reference, Val(S), frequency, wavelength)
    lambda_reference = _fits_parent_to_x(parent_reference, Val(X), Val(P), frequency, wavelength)
    ruling, order, incidence, refractive_index, refractive_derivative, grating_tilt, detector_tilt = parameters
    alpha = deg2rad(incidence)
    epsilon = deg2rad(grating_tilt)
    theta = deg2rad(detector_tilt)
    grating = ruling * order / cos(epsilon)
    beta_argument = grating * lambda_reference - refractive_index * sin(alpha)
    abs(beta_argument) <= 1 || throw(DomainError(beta_argument, "FITS grism reference wavelength lies outside the inverse-sine domain"))
    beta_reference = asin(beta_argument)
    denominator = grating - refractive_derivative * sin(alpha)
    iszero(denominator) && throw(ArgumentError("FITS grism parameters define a singular wavelength equation"))
    constant = (refractive_index - refractive_derivative * lambda_reference) * sin(alpha)
    gamma_reference = -tan(theta)
    lambda_derivative = ForwardDiff.derivative(reference) do value
        parent = _fits_spectral_to_parent(value, Val(S), frequency, wavelength)
        _fits_parent_to_x(parent, Val(X), Val(P), frequency, wavelength)
    end
    linear_scale = lambda_derivative * denominator / (cos(beta_reference) * cos(theta)^2)
    all(isfinite, (lambda_reference, gamma_reference, beta_reference, beta_reference + theta, constant, denominator, linear_scale, scale)) ||
        throw(ArgumentError("FITS grism parameters do not define a finite transform"))
    iszero(linear_scale) && throw(ArgumentError("FITS grism parameters define a singular reference derivative"))
    return FITSGrismTransform{X,P,S,A,typeof(reference)}(gamma_reference, beta_reference, beta_reference + theta, constant,
        denominator, frequency, wavelength, scale, linear_scale)
end

@inline function (transform::FITSGrismTransform{X,P,S,A})(coords::StaticVector{1}) where {X,P,S,A}
    offset = only(coords)
    gamma = oftype(offset, transform.gamma_reference) + offset
    beta = atan(gamma) + oftype(offset, transform.beta_offset)
    wavelength = (sin(beta) + oftype(offset, transform.constant)) / oftype(offset, transform.denominator)
    parent = _fits_x_to_parent(wavelength, Val(X), Val(P), oftype(offset, transform.rest_frequency), oftype(offset, transform.rest_wavelength))
    world = _fits_parent_to_spectral(parent, Val(S), oftype(offset, transform.rest_frequency), oftype(offset, transform.rest_wavelength))
    return SVector(world / oftype(world, transform.output_scale))
end

struct FITSGrismInverse{X,P,S,A,F} <: AbstractTransform{1,1}
    forward::F
end

@inline function (transform::FITSGrismInverse{X,P,S,A})(coords::StaticVector{1}) where {X,P,S,A}
    forward = transform.forward
    world = only(coords) * oftype(only(coords), forward.output_scale)
    parent = _fits_spectral_to_parent(world, Val(S), oftype(world, forward.rest_frequency), oftype(world, forward.rest_wavelength))
    wavelength = _fits_parent_to_x(parent, Val(X), Val(P), oftype(world, forward.rest_frequency), oftype(world, forward.rest_wavelength))
    sine_beta = wavelength * oftype(world, forward.denominator) - oftype(world, forward.constant)
    abs(sine_beta) <= one(sine_beta) || return SVector(oftype(sine_beta, NaN))
    gamma = tan(asin(sine_beta) - oftype(world, forward.beta_offset))
    return SVector(gamma - oftype(world, forward.gamma_reference))
end

inverse_transform(transform::FITSGrismTransform{X,P,S,A}) where {X,P,S,A} = FITSGrismInverse{X,P,S,A,typeof(transform)}(transform)
inverse_transform(transform::FITSGrismInverse) = transform.forward
