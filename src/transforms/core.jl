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
transform_name(::AbstractTransform) = nothing
input_names(::AbstractTransform) = nothing
output_names(::AbstractTransform) = nothing

function _dense_dependencies(::Val{Nout}, ::Val{Nin}) where {Nin,Nout}
    return SMatrix{Nout,Nin,Bool}(ntuple(_ -> true, Nin * Nout))
end

"""
    transform_dependencies(transform)

Return a boolean output-by-input matrix describing structural coordinate
dependencies. External transform packages may specialize this method; the
default conservatively assumes every output depends on every input.
"""
transform_dependencies(::AbstractTransform{Nin,Nout}) where {Nin,Nout} = _dense_dependencies(Val(Nout), Val(Nin))

function _dependency_matrix(dependencies, ::Val{Nout}, ::Val{Nin}, context) where {Nin,Nout}
    dependencies isa AbstractMatrix || throw(ArgumentError("$context returned $(typeof(dependencies)), expected an AbstractMatrix"))
    size(dependencies) == (Nout, Nin) || throw(DimensionMismatch("$context has size $(size(dependencies)), expected ($Nout, $Nin)"))
    try
        return SMatrix{Nout,Nin,Bool}(dependencies)
    catch error
        error isa InexactError || error isa MethodError || rethrow()
        throw(ArgumentError("$context entries must be boolean-convertible"))
    end
end

function _transform_dependencies_checked(transform::AbstractTransform{Nin,Nout}) where {Nin,Nout}
    return _dependency_matrix(transform_dependencies(transform), Val(Nout), Val(Nin), "transform_dependencies(::$(typeof(transform)))")
end

axis_correlation_matrix(transform::AbstractTransform) = _transform_dependencies_checked(transform)

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

"""
    transform_jacobian(transform, coords)

Return the Jacobian of a transform with respect to its inputs. External
transform packages may extend this method to provide an analytical Jacobian.
"""
function _forwarddiff_jacobian(transform::AbstractTransform{Nin}, coords::StaticVector{Nin}) where {Nin}
    T = _evaluation_type(eltype(coords))
    return ForwardDiff.jacobian(transform, SVector{Nin,T}(coords))
end

transform_jacobian(transform::AbstractTransform, coords::StaticVector) = _forwarddiff_jacobian(transform, coords)

function _transform_jacobian_checked(transform::AbstractTransform{Nin,Nout}, coords::StaticVector{Nin}) where {Nin,Nout}
    jacobian = transform_jacobian(transform, coords)
    jacobian isa AbstractMatrix || throw(ArgumentError("transform_jacobian(::$(typeof(transform)), ::$(typeof(coords))) returned $(typeof(jacobian)), expected an AbstractMatrix"))
    size(jacobian) == (Nout, Nin) || throw(DimensionMismatch("Jacobian of $(typeof(transform)) has size $(size(jacobian)), expected ($Nout, $Nin)"))
    return SMatrix{Nout,Nin}(jacobian)
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
transform_jacobian(transform::WithInverse, coords::StaticVector) = _transform_jacobian_checked(transform.forward, coords)
transform_name(transform::WithInverse) = transform_name(transform.forward)
input_names(transform::WithInverse) = input_names(transform.forward)
output_names(transform::WithInverse) = output_names(transform.forward)
