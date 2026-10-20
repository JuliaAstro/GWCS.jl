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

One-dimensional polynomial with coefficients in increasing power order.
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

Triangular two-dimensional polynomial. Matrix element `[i + 1, j + 1]` is the
coefficient of `x^i * y^j`; entries where `i + j` exceeds the degree are ignored.
"""
struct Polynomial2D{N,T,L,D,W} <: AbstractTransform{2,1}
    coefficients::SMatrix{N,N,T,L}
    domain::D
    window::W
end

function Polynomial2D(coefficients::AbstractMatrix; domain=nothing, window=nothing)
    size(coefficients, 1) == size(coefficients, 2) || throw(ArgumentError("coefficients must be a square matrix"))
    N = size(coefficients, 1)
    N > 0 || throw(ArgumentError("polynomial requires at least one coefficient"))
    (domain === nothing) == (window === nothing) || throw(ArgumentError("domain and window must be provided together"))
    T = eltype(coefficients)
    return Polynomial2D{N,T,N * N,typeof(domain),typeof(window)}(SMatrix{N,N,T}(coefficients), domain, window)
end

degree(::Polynomial2D{N}) where {N} = N - 1

@inline function (transform::Polynomial2D{N})(coords::StaticVector{2}) where {N}
    T = _evaluation_type(eltype(coords))
    x = _map_domain(T(coords[1]), transform.domain === nothing ? nothing : transform.domain[1], transform.window === nothing ? nothing : transform.window[1])
    y = _map_domain(T(coords[2]), transform.domain === nothing ? nothing : transform.domain[2], transform.window === nothing ? nothing : transform.window[2])
    result = zero(x + y + T(transform.coefficients[1, 1]))
    @inbounds for j in (N - 1):-1:0
        max_i = N - 1 - j
        xpoly = T(transform.coefficients[max_i + 1, j + 1])
        for i in (max_i - 1):-1:0
            xpoly = muladd(xpoly, x, T(transform.coefficients[i + 1, j + 1]))
        end
        result = muladd(result, y, xpoly)
    end
    return SVector(result)
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
function _linear_form(transform::Polynomial2D{N}) where {N}
    c = transform.coefficients
    cx, cy = N > 1 ? (c[2, 1], c[1, 2]) : (zero(c[1, 1]), zero(c[1, 1]))
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

"""A TPV plane polynomial using FITS sequent-polynomial term ordering."""
struct TPVTransform{NX,NY,T,TX,TY} <: AbstractTransform{2,2}
    xcoeff::NTuple{NX,T}
    ycoeff::NTuple{NY,T}
    xterms::TX
    yterms::TY
end

function TPVTransform(xcoeff, ycoeff)
    xvalues, yvalues = Tuple(float.(xcoeff)), Tuple(float.(ycoeff))
    length(xvalues) <= 40 && length(yvalues) <= 40 || throw(ArgumentError("TPV supports coefficient indices 0 through 39"))
    T = isempty(xvalues) && isempty(yvalues) ? Float64 : promote_type(map(typeof, (xvalues..., yvalues...))...)
    xterms = ntuple(index -> _tpd_term_descriptor(index - 1), length(xvalues))
    yterms = ntuple(index -> _tpd_term_descriptor(index - 1), length(yvalues))
    return TPVTransform{length(xvalues),length(yvalues),T,typeof(xterms),typeof(yterms)}(T.(xvalues), T.(yvalues), xterms, yterms)
end

@inline function _tpd_evaluate(coefficients, terms, x, y)
    radius = hypot(x, y)
    result = zero(x + y)
    for index in eachindex(coefficients)
        radial, xpower, ypower = terms[index]
        value = radial ? radius^xpower : x^xpower * y^ypower
        result = muladd(typeof(result)(coefficients[index]), value, result)
    end
    return result
end

@inline function (transform::TPVTransform)(coords::StaticVector{2})
    x, y = coords
    return SVector(_tpd_evaluate(transform.xcoeff, transform.xterms, x, y), _tpd_evaluate(transform.ycoeff, transform.yterms, y, x))
end

struct _FITSTPDAxis{M,O,S,A,T,R}
    axes::NTuple{M,Int}
    offset::O
    scale::S
    auxiliary::A
    terms::T
    reverse_terms::R
    docorr::Bool
end

function _FITSTPDAxis(axes::NTuple{M,Int}, offset, scale, auxiliary, terms, reverse_terms, docorr) where {M}
    return _FITSTPDAxis{M,typeof(offset),typeof(scale),typeof(auxiliary),typeof(terms),typeof(reverse_terms)}(axes, offset, scale, auxiliary, terms, reverse_terms, docorr)
end

_FITSTPDAxis(axes, offset, scale, auxiliary, terms, docorr) = _FITSTPDAxis(axes, offset, scale, auxiliary, terms, (), docorr)

struct _FITSTPDTransform{N,M<:Tuple} <: AbstractTransform{N,N}
    models::M
end

function _FITSTPDTransform(models::NTuple{N,Any}) where {N}
    return _FITSTPDTransform{N,typeof(models)}(models)
end

@inline function _fits_tpd_value(model::_FITSTPDAxis{M}, coords::StaticVector) where {M}
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
    radius = hypot(u, v)
    result = zero(T)
    for (coefficient, radial, xpower, ypower) in model.terms
        value = radial ? radius^xpower : u^xpower * v^ypower
        result = muladd(T(coefficient), value, result)
    end
    return result
end

@inline _fits_tpd_output(::Nothing, coords, axis) = coords[axis]

@inline function _fits_tpd_output(model::_FITSTPDAxis, coords, axis)
    value = _fits_tpd_value(model, coords)
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
    isempty(model.reverse_terms) && return (coords[axis], false)
    reverse_model = _FITSTPDAxis(model.axes, model.offset, model.scale, model.auxiliary, model.reverse_terms, (), model.docorr)
    return (_fits_tpd_output(reverse_model, coords, axis), true)
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
