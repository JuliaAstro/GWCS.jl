@inline function _map_domain(value, domain::Nothing, window)
    return value
end

@inline function _map_domain(value, domain::Tuple, window::Tuple)
    T = typeof(value)
    domain_low, domain_high = T(domain[1]), T(domain[2])
    window_low, window_high = T(window[1]), T(window[2])
    domain_mid = (domain_low + domain_high) / T(2)
    window_mid = (window_low + window_high) / T(2)
    return (value - domain_mid) * (window_high - window_low) / (domain_high - domain_low) + window_mid
end

"""
    Polynomial1D(coefficients; domain=nothing, window=nothing)

The one-dimensional polynomial with `coefficients` in increasing powers, `c[1] + c[2] * x + c[3] * x^2 + ...`. With
`domain` and `window`, two `(low, high)` pairs given together, the input is first mapped linearly from `domain` to
`window`. It has no inverse.

Python equivalent: `astropy.modeling.models.Polynomial1D`.
"""
struct Polynomial1D{N,T,D,W} <: AbstractTransform{1,1}
    coefficients::SVector{N,T}
    domain::D
    window::W
end

function Polynomial1D(coefficients::AbstractVector; domain=nothing, window=nothing)
    N = length(coefficients)
    N > 0 || throw(ArgumentError("polynomial requires at least one coefficient"))
    (domain === nothing) == (window === nothing) || throw(ArgumentError("domain and window must be provided together"))
    T = eltype(coefficients)
    return Polynomial1D{N,T,typeof(domain),typeof(window)}(SVector{N,T}(coefficients), domain, window)
end

"""
    degree(polynomial::Union{Polynomial1D,Polynomial2D})

Return the degree a polynomial's coefficients are stored for, whether or not its highest coefficients are zero. It is
a compile-time constant for static coefficient storage.

Python equivalent: the `degree` attribute of `astropy.modeling.models.Polynomial1D` and `Polynomial2D`.
"""
degree(::Polynomial1D{N}) where {N} = N - 1

@inline function (transform::Polynomial1D{N})(coords::StaticVector{1}) where {N}
    T = _evaluation_type(eltype(coords))
    x = _map_domain(T(coords[1]), transform.domain, transform.window)
    value = T(transform.coefficients[N])
    @inbounds for i in (N - 1):-1:1
        value = muladd(value, x, T(transform.coefficients[i]))
    end
    return SVector(value)
end

"""
    Polynomial2D(coefficients; domain=nothing, window=nothing)
    Polynomial2D{T,C,D,W}(coefficients, domain, window)

The two-dimensional polynomial of total [`degree`](@ref) `n`, given by a square `(n + 1) × (n + 1)` coefficient matrix
whose element `[i + 1, j + 1]` multiplies `x^i * y^j`; elements with `i + j > n` are ignored. `domain` and `window`,
each a pair of `(low, high)` intervals for `x` and `y`, map the inputs as for [`Polynomial1D`](@ref). It has no
inverse.

The first form stores the coefficients as an `SMatrix`, so the degree is part of the type and evaluation unrolls. The
parametric form chooses the coefficient storage `C`, with `D` and `W` the types of `domain` and `window`. With
`C = Matrix{T}`, polynomials of every degree share one type and one compiled kernel, at two to three times the
evaluation cost.

Python equivalent: `astropy.modeling.models.Polynomial2D`, whose coefficients are named `c{i}_{j}` for the same terms.
"""
struct Polynomial2D{T,C<:AbstractMatrix{T},D,W} <: AbstractTransform{2,1}
    coefficients::C
    domain::D
    window::W
end

function Polynomial2D(coefficients::AbstractMatrix; domain=nothing, window=nothing)
    size(coefficients, 1) == size(coefficients, 2) || throw(ArgumentError("coefficients must be a square matrix"))
    N = size(coefficients, 1)
    N > 0 || throw(ArgumentError("polynomial requires at least one coefficient"))
    (domain === nothing) == (window === nothing) || throw(ArgumentError("domain and window must be provided together"))
    T = eltype(coefficients)
    return Polynomial2D{T,SMatrix{N,N,T,N * N},typeof(domain),typeof(window)}(coefficients, domain, window)
end

degree(transform::Polynomial2D) = size(transform.coefficients, 1) - 1

@inline function (transform::Polynomial2D)(coords::StaticVector{2})
    T = _evaluation_type(eltype(coords))
    x = _map_domain(T(coords[1]), transform.domain === nothing ? nothing : transform.domain[1], transform.window === nothing ? nothing : transform.window[1])
    y = _map_domain(T(coords[2]), transform.domain === nothing ? nothing : transform.domain[2], transform.window === nothing ? nothing : transform.window[2])
    return SVector(_polynomial2d_value(transform.coefficients, x, y))
end

# Horner's rule over a triangular coefficient matrix whose element `[i + 1, j + 1]` multiplies `x^i * y^j`. For an
# SMatrix the size is a constant, so the loops unroll.
@inline function _polynomial2d_value(coefficients::AbstractMatrix, x::T, y::T) where {T}
    N = size(coefficients, 1)
    result = zero(x + y + T(coefficients[1, 1]))
    @inbounds for j in (N - 1):-1:0
        max_i = N - 1 - j
        xpoly = T(coefficients[max_i + 1, j + 1])
        for i in (max_i - 1):-1:0
            xpoly = muladd(xpoly, x, T(coefficients[i + 1, j + 1]))
        end
        result = muladd(result, y, xpoly)
    end
    return result
end

# Python GWCS's `Mapping((0, 1, 0, 1)) | Polynomial2D & Polynomial2D` idiom as one leaf, which `compose` builds.
struct _Polynomial2DPair{X<:Polynomial2D,Y<:Polynomial2D} <: AbstractTransform{2,2}
    x::X
    y::Y
end

@inline (transform::_Polynomial2DPair)(coords::StaticVector{2}) = SVector(only(transform.x(coords)), only(transform.y(coords)))

# astropy's `poly_map_domain` written as x -> scale * x + offset.
_domain_affine(domain::Nothing, window) = (1, 0)
_domain_affine(domain, window) = ((window[2] - window[1]) / (domain[2] - domain[1]), (window[1] * domain[2] - window[2] * domain[1]) / (domain[2] - domain[1]))

# The gradient row and constant of a polynomial of degree at most one.
function _linear_form(transform::Polynomial2D)
    c = transform.coefficients
    cx, cy = size(c, 1) > 1 ? (c[2, 1], c[1, 2]) : (zero(c[1, 1]), zero(c[1, 1]))
    sx, ox = _domain_affine(transform.domain === nothing ? nothing : transform.domain[1], transform.window === nothing ? nothing : transform.window[1])
    sy, oy = _domain_affine(transform.domain === nothing ? nothing : transform.domain[2], transform.window === nothing ? nothing : transform.window[2])
    return SVector(cx * sx, cy * sy), c[1, 1] + cx * ox + cy * oy
end

"""A two-dimensional SIP correction `(u, v) -> (u + A(u,v), v + B(u,v))`."""
struct SIPTransform{A<:Polynomial2D,B<:Polynomial2D} <: AbstractTransform{2,2}
    a::A
    b::B
end

function _tpd_term_descriptor(index)
    remaining = index
    degree = 0
    while true
        monomials = degree + 1
        if remaining < monomials
            return (false, degree - remaining, remaining)
        end
        remaining -= monomials
        if isodd(degree)
            iszero(remaining) && return (true, degree, 0)
            remaining -= 1
        end
        degree += 1
    end
end

# Sequent TPD coefficients (FITS order, index 0 first) as `(coefficient, radial, xpower, ypower)` terms, skipping zeros.
function _fits_tpd_compiled_terms(coefficients, ::Val{M}; radial=true) where {M}
    terms = Tuple{Float64,Bool,Int,Int}[]
    for (position, coefficient) in enumerate(coefficients)
        isradial, xpower, ypower = _tpd_term_descriptor(position - 1)
        M == 1 && (isradial || !iszero(ypower)) && continue
        (!radial && isradial) || iszero(coefficient) || push!(terms, (Float64(coefficient), isradial, xpower, ypower))
    end
    return Tuple(terms)
end

# Like WCSLIB's tpd1-tpd9, evaluate every coefficient up to the degree by Horner's rule rather than each term by
# integer powers: `coefficients[i + 1, j + 1]` multiplies u^i v^j and `radial[k]` multiplies r^(2k - 1).
function _tpd_polynomial(terms)
    degree = maximum((xpower + ypower for (_, radial, xpower, ypower) in terms if !radial); init=0)
    radial_count = maximum(((xpower + 1) ÷ 2 for (_, radial, xpower, _) in terms if radial); init=0)
    coefficients = zeros(degree + 1, degree + 1)
    radial = zeros(radial_count)
    for (coefficient, isradial, xpower, ypower) in terms
        isradial ? (radial[(xpower + 1) ÷ 2] += coefficient) : (coefficients[xpower + 1, ypower + 1] += coefficient)
    end
    return (SMatrix{degree + 1,degree + 1}(coefficients), Tuple(radial))
end

@inline function _tpd_polynomial_value(polynomial, u, v)
    coefficients, radial = polynomial
    result = _polynomial2d_value(coefficients, u, v)
    isempty(radial) && return result
    s = muladd(u, u, v * v)
    return muladd(sqrt(s), evalpoly(s, radial), result)
end

"""
    TPVTransform(xcoeff, ycoeff)

TPV plane polynomial in FITS sequent order, including radial terms: `xcoeff` holds `PV1_0`, `PV1_1`, ... and gives
ξ(x, y); `ycoeff` holds `PV2_*` and gives η(y, x).
"""
struct TPVTransform{X,Y} <: AbstractTransform{2,2}
    x::X
    y::Y
    # An inner constructor suppresses the default outer one, which would take the place of the coefficient constructor.
    TPVTransform{X,Y}(x, y) where {X,Y} = new{X,Y}(x, y)
end

function TPVTransform(xcoeff, ycoeff)
    length(xcoeff) <= 40 && length(ycoeff) <= 40 || throw(ArgumentError("TPV supports coefficient indices 0 through 39"))
    x, y = _tpd_polynomial(_fits_tpd_compiled_terms(xcoeff, Val(2))), _tpd_polynomial(_fits_tpd_compiled_terms(ycoeff, Val(2)))
    return TPVTransform{typeof(x),typeof(y)}(x, y)
end

@inline function (transform::TPVTransform)(coords::StaticVector{2})
    T = _evaluation_type(eltype(coords))
    x, y = T(coords[1]), T(coords[2])
    return SVector(_tpd_polynomial_value(transform.x, x, y), _tpd_polynomial_value(transform.y, y, x))
end

struct _FITSTPDAxis{M,O,S,A,P,R}
    axes::NTuple{M,Int}
    offset::O
    scale::S
    auxiliary::A
    polynomial::P
    reverse::R
    docorr::Bool
    # An inner constructor suppresses the default outer one, which would bypass the term conversion below.
    _FITSTPDAxis{M,O,S,A,P,R}(axes, offset, scale, auxiliary, polynomial, reverse, docorr) where {M,O,S,A,P,R} = new{M,O,S,A,P,R}(axes, offset, scale, auxiliary, polynomial, reverse, docorr)
end

# Terms are `(coefficient, radial, xpower, ypower)` tuples; the reverse polynomial is `nothing` when there are none.
function _FITSTPDAxis(axes::NTuple{M,Int}, offset, scale, auxiliary, terms, reverse_terms, docorr) where {M}
    polynomial = _tpd_polynomial(terms)
    reverse = isempty(reverse_terms) ? nothing : _tpd_polynomial(reverse_terms)
    return _FITSTPDAxis{M,typeof(offset),typeof(scale),typeof(auxiliary),typeof(polynomial),typeof(reverse)}(axes, offset, scale, auxiliary, polynomial, reverse, docorr)
end

_FITSTPDAxis(axes, offset, scale, auxiliary, terms, docorr) = _FITSTPDAxis(axes, offset, scale, auxiliary, terms, (), docorr)

struct _FITSTPDTransform{N,M<:Tuple} <: AbstractTransform{N,N}
    models::M
end

function _FITSTPDTransform(models::NTuple{N,Any}) where {N}
    return _FITSTPDTransform{N,typeof(models)}(models)
end

@inline function _fits_tpd_value(model::_FITSTPDAxis{M}, polynomial, coords::StaticVector) where {M}
    T = _evaluation_type(eltype(coords))
    variables = SVector{M,T}(ntuple(Val(M)) do index
        (T(coords[model.axes[index]]) - T(model.offset[index])) * T(model.scale[index])
    end)
    u = variables[1]
    v = M == 1 ? zero(T) : variables[2]
    if model.auxiliary !== nothing
        transformed = SMatrix{2,3,T}(model.auxiliary) * SVector(one(T), u, v)
        u, v = transformed
    end
    return _tpd_polynomial_value(polynomial, u, v)
end

@inline _fits_tpd_output(::Nothing, coords, axis) = coords[axis]

@inline function _fits_tpd_output(model::_FITSTPDAxis, coords, axis)
    value = _fits_tpd_value(model, model.polynomial, coords)
    return model.docorr ? coords[axis] + value : value
end

@inline _fits_tpd_outputs(::Tuple{}, coords, axis, ::Type{T}) where {T} = ()

@inline function _fits_tpd_outputs(models::Tuple, coords, axis, ::Type{T}) where {T}
    output = T(_fits_tpd_output(first(models), coords, axis))
    return (output, _fits_tpd_outputs(Base.tail(models), coords, axis + 1, T)...)
end

@inline function (transform::_FITSTPDTransform{N})(coords::StaticVector{N}) where {N}
    T = _evaluation_type(eltype(coords))
    values = SVector{N,T}(coords)
    return SVector{N,T}(_fits_tpd_outputs(transform.models, values, 1, T))
end

@inline _fits_tpd_reverse_output(::Nothing, coords, axis) = (coords[axis], false)

@inline function _fits_tpd_reverse_output(model::_FITSTPDAxis, coords, axis)
    model.reverse === nothing && return (coords[axis], false)
    value = _fits_tpd_value(model, model.reverse, coords)
    return (model.docorr ? coords[axis] + value : value, true)
end

@inline _fits_tpd_reverse_outputs(::Tuple{}, coords, axis, ::Type{T}) where {T} = ((), false)

@inline function _fits_tpd_reverse_outputs(models::Tuple, coords, axis, ::Type{T}) where {T}
    output, used = _fits_tpd_reverse_output(first(models), coords, axis)
    tail, tail_used = _fits_tpd_reverse_outputs(Base.tail(models), coords, axis + 1, T)
    return ((T(output), tail...), used || tail_used)
end

@inline function _inverse_initialize(transform::_FITSTPDTransform{N}, coords::StaticVector{N}) where {N}
    T = _evaluation_type(eltype(coords))
    values = SVector{N,T}(coords)
    outputs, used = _fits_tpd_reverse_outputs(transform.models, values, 1, T)
    return (SVector{N,T}(outputs), used)
end

@inline function (transform::SIPTransform)(coords::StaticVector{2})
    return coords + SVector(only(transform.a(coords)), only(transform.b(coords)))
end
