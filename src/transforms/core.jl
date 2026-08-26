"""
    AbstractTransform{Nin,Nout}

Abstract supertype for transformations with statically known input and output
dimensions.
"""
abstract type AbstractTransform{Nin,Nout} end

"""Marker returned when a transform has no compiled inverse."""
struct NoInverse end

input_dimension(::AbstractTransform{Nin,Nout}) where {Nin,Nout} = Nin
output_dimension(::AbstractTransform{Nin,Nout}) where {Nin,Nout} = Nout

"""Return the analytical or explicitly supplied inverse, or `NoInverse()`."""
inverse_transform(::AbstractTransform) = NoInverse()

has_inverse(transform::AbstractTransform) = !(inverse_transform(transform) isa NoInverse)

"""Compile a transform into its execution representation."""
compile_transform(transform::AbstractTransform) = transform

_evaluation_type(::Type{T}) where {T<:Integer} = Float64
_evaluation_type(::Type{T}) where {T<:Rational} = Float64
_evaluation_type(::Type{T}) where {T} = T

"""
    WithInverse(forward, backward)

Pair a transform with an explicit inverse supplied by metadata or a user.
"""
struct WithInverse{Nin,Nout,F<:AbstractTransform{Nin,Nout},I<:AbstractTransform{Nout,Nin}} <: AbstractTransform{Nin,Nout}
    forward::F
    backward::I
end

@inline (transform::WithInverse)(coords::StaticVector) = transform.forward(coords)

inverse_transform(transform::WithInverse) = WithInverse(transform.backward, transform.forward)
compile_transform(transform::WithInverse) = WithInverse(compile_transform(transform.forward), compile_transform(transform.backward))
