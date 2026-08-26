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
