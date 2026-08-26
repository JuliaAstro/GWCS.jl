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

"""
    inverse_transform(transform)

Return the analytical or explicitly supplied inverse, or `NoInverse()`.
External transform packages may extend this method for native
`AbstractTransform` subtypes.
"""
inverse_transform(::AbstractTransform) = NoInverse()

function _inverse_transform_checked(transform::AbstractTransform)
    inverse = inverse_transform(transform)
    inverse isa NoInverse && return inverse
    inverse isa AbstractTransform || throw(ArgumentError("inverse_transform(::$(typeof(transform))) returned $(typeof(inverse)), expected AbstractTransform or NoInverse"))
    input_dimension(inverse) == output_dimension(transform) && output_dimension(inverse) == input_dimension(transform) || throw(DimensionMismatch("inverse of $(typeof(transform)) has dimensions $(input_dimension(inverse)) -> $(output_dimension(inverse)), expected $(output_dimension(transform)) -> $(input_dimension(transform))"))
    return inverse
end

has_inverse(transform::AbstractTransform) = !(_inverse_transform_checked(transform) isa NoInverse)

"""
    compile_transform(transform)

Compile a transform into a native execution representation with identical
dimensions. External transform packages may extend this method to replace a
construction-oriented transform with an optimized `AbstractTransform`.
"""
compile_transform(transform::AbstractTransform) = transform

function _compile_transform_checked(transform::AbstractTransform)
    compiled = compile_transform(transform)
    compiled isa AbstractTransform || throw(ArgumentError("compile_transform(::$(typeof(transform))) returned $(typeof(compiled)), expected AbstractTransform"))
    input_dimension(compiled) == input_dimension(transform) && output_dimension(compiled) == output_dimension(transform) || throw(DimensionMismatch("compiled form of $(typeof(transform)) has dimensions $(input_dimension(compiled)) -> $(output_dimension(compiled)), expected $(input_dimension(transform)) -> $(output_dimension(transform))"))
    return compiled
end

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
compile_transform(transform::WithInverse) = WithInverse(_compile_transform_checked(transform.forward), _compile_transform_checked(transform.backward))
