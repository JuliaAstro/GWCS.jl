"""
    CallableTransform{Nin,Nout}(forward; inverse=nothing, jacobian=nothing, dependencies=nothing, inverse_dependencies=nothing)
    CallableTransform(forward, Nin, Nout; kwargs...)

A transform that calls `forward(coords...)` with `Nin` numbers. The callable returns a tuple or vector of `Nout`
values, or a number when `Nout` is 1.

`inverse` is a callable of the same form for the inverse. `jacobian(coords...)` returns the `Nout × Nin` Jacobian;
without it, ForwardDiff.jl computes the Jacobian. `dependencies` and `inverse_dependencies`, `Bool` matrices like
[`transform_dependencies`](@ref)'s, declare which outputs depend on which inputs.

A `DomainError` thrown by a callable gives `NaN`, as GWCS.jl transforms do where a coordinate has no image. It cannot
be written to metadata.
"""
struct CallableTransform{Nin,Nout,F,I,J,D,ID} <: AbstractTransform{Nin,Nout}
    forward::F
    inverse::I
    jacobian::J
    dependencies::D
    inverse_dependencies::ID
end

struct InverseCallableTransform{Nin,Nout,T<:CallableTransform{Nout,Nin}} <: AbstractTransform{Nin,Nout}
    transform::T
end

function _callable_dependencies(dependencies, ::Val{Nin}, ::Val{Nout}, name) where {Nin,Nout}
    dependencies === nothing && return nothing
    return _dependency_matrix(dependencies, Val(Nout), Val(Nin), name)
end

function CallableTransform{Nin,Nout}(forward; inverse=nothing, jacobian=nothing, dependencies=nothing, inverse_dependencies=nothing) where {Nin,Nout}
    Nin >= 0 || throw(ArgumentError("input dimension must be non-negative"))
    Nout > 0 || throw(ArgumentError("output dimension must be positive"))
    inverse === nothing && inverse_dependencies !== nothing && throw(ArgumentError("inverse_dependencies requires an inverse callable"))
    forward_dependencies = _callable_dependencies(dependencies, Val(Nin), Val(Nout), "callable dependencies")
    backward_dependencies = _callable_dependencies(inverse_dependencies, Val(Nout), Val(Nin), "callable inverse dependencies")
    return CallableTransform{Nin,Nout,typeof(forward),typeof(inverse),typeof(jacobian),typeof(forward_dependencies),typeof(backward_dependencies)}(forward, inverse, jacobian, forward_dependencies, backward_dependencies)
end

function CallableTransform(forward, nin::Integer, nout::Integer; inverse=nothing, jacobian=nothing, dependencies=nothing, inverse_dependencies=nothing)
    return CallableTransform{Int(nin),Int(nout)}(forward; inverse, jacobian, dependencies, inverse_dependencies)
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

# User functions may throw DomainError (for example Base `sqrt` of a negative number) where GWCS
# transforms return NaN, so convert it; numerical-inverse steps routinely leave a function's domain.
@inline function _callable_evaluate(f, ::Val{N}, coords) where {N}
    try
        return _callable_output(Val(N), f(coords...))
    catch error
        error isa DomainError || rethrow()
        return SVector{N}(ntuple(_ -> _evaluation_type(eltype(coords))(NaN), N))
    end
end

@inline (transform::CallableTransform{Nin,Nout})(coords::StaticVector{Nin}) where {Nin,Nout} = _callable_evaluate(transform.forward, Val(Nout), coords)

inverse_transform(::CallableTransform{Nin,Nout,F,Nothing,J,D,ID}) where {Nin,Nout,F,J,D,ID} = NoInverse()

function inverse_transform(transform::CallableTransform{Nin,Nout}) where {Nin,Nout}
    return InverseCallableTransform{Nout,Nin,typeof(transform)}(transform)
end

@inline (transform::InverseCallableTransform{Nin,Nout})(coords::StaticVector{Nin}) where {Nin,Nout} = _callable_evaluate(transform.transform.inverse, Val(Nout), coords)

inverse_transform(transform::InverseCallableTransform) = transform.transform

function _callable_jacobian(::Val{1}, ::Val{1}, result::Number)
    return SMatrix{1,1}(result)
end

function _callable_jacobian(::Val{Nin}, ::Val{Nout}, result) where {Nin,Nout}
    result isa AbstractMatrix || throw(ArgumentError("callable Jacobian returned $(typeof(result)), expected an AbstractMatrix"))
    size(result) == (Nout, Nin) || throw(DimensionMismatch("callable Jacobian has size $(size(result)), expected ($Nout, $Nin)"))
    return SMatrix{Nout,Nin}(result)
end

function transform_jacobian(transform::CallableTransform{Nin,Nout,F,I,Nothing,D,ID}, coords::StaticVector{Nin}) where {Nin,Nout,F,I,D,ID}
    return _forwarddiff_jacobian(transform, coords)
end

function transform_jacobian(transform::CallableTransform{Nin,Nout}, coords::StaticVector{Nin}) where {Nin,Nout}
    result = transform.jacobian(coords...)
    return _callable_jacobian(Val(Nin), Val(Nout), result)
end
