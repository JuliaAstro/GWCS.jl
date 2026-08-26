"""
    CallableTransform{Nin,Nout}(forward; inverse=nothing, jacobian=nothing)
    CallableTransform(forward, nin, nout; inverse=nothing, jacobian=nothing)

Wrap positional scalar callables as a statically dimensioned transform.
One-output callables may return a scalar; otherwise callables return a tuple or
vector with the declared number of outputs.
"""
struct CallableTransform{Nin,Nout,F,I,J} <: AbstractTransform{Nin,Nout}
    forward::F
    inverse::I
    jacobian::J
end

struct InverseCallableTransform{Nin,Nout,T<:CallableTransform{Nout,Nin}} <: AbstractTransform{Nin,Nout}
    transform::T
end

function CallableTransform{Nin,Nout}(forward; inverse=nothing, jacobian=nothing) where {Nin,Nout}
    Nin >= 0 || throw(ArgumentError("input dimension must be non-negative"))
    Nout > 0 || throw(ArgumentError("output dimension must be positive"))
    return CallableTransform{Nin,Nout,typeof(forward),typeof(inverse),typeof(jacobian)}(forward, inverse, jacobian)
end

function CallableTransform(forward, nin::Integer, nout::Integer; inverse=nothing, jacobian=nothing)
    return CallableTransform{Int(nin),Int(nout)}(forward; inverse, jacobian)
end

@inline function _callable_output(::Val{N}, result) where {N}
    if result isa Number
        N == 1 || throw(DimensionMismatch("callable transform returned one scalar, expected $N outputs"))
        return SVector(result)
    end
    result isa Tuple || result isa AbstractVector || throw(ArgumentError("callable transform returned $(typeof(result)), expected a scalar, tuple, or AbstractVector"))
    length(result) == N || throw(DimensionMismatch("callable transform returned $(length(result)) outputs, expected $N"))
    return SVector{N}(result)
end

@inline function (transform::CallableTransform{Nin,Nout})(coords::StaticVector{Nin}) where {Nin,Nout}
    return _callable_output(Val(Nout), transform.forward(coords...))
end

inverse_transform(::CallableTransform{Nin,Nout,F,Nothing,J}) where {Nin,Nout,F,J} = NoInverse()

function inverse_transform(transform::CallableTransform{Nin,Nout}) where {Nin,Nout}
    return InverseCallableTransform{Nout,Nin,typeof(transform)}(transform)
end

@inline function (transform::InverseCallableTransform{Nin,Nout})(coords::StaticVector{Nin}) where {Nin,Nout}
    return _callable_output(Val(Nout), transform.transform.inverse(coords...))
end

inverse_transform(transform::InverseCallableTransform) = transform.transform

function _callable_jacobian(::Val{1}, ::Val{1}, result::Number)
    return SMatrix{1,1}(result)
end

function _callable_jacobian(::Val{Nin}, ::Val{Nout}, result) where {Nin,Nout}
    result isa AbstractMatrix || throw(ArgumentError("callable Jacobian returned $(typeof(result)), expected an AbstractMatrix"))
    size(result) == (Nout, Nin) || throw(DimensionMismatch("callable Jacobian has size $(size(result)), expected ($Nout, $Nin)"))
    return SMatrix{Nout,Nin}(result)
end

function transform_jacobian(transform::CallableTransform{Nin,Nout,F,I,Nothing}, coords::StaticVector{Nin}) where {Nin,Nout,F,I}
    return _forwarddiff_jacobian(transform, coords)
end

function transform_jacobian(transform::CallableTransform{Nin,Nout}, coords::StaticVector{Nin}) where {Nin,Nout}
    result = transform.jacobian(coords...)
    return _callable_jacobian(Val(Nin), Val(Nout), result)
end
