# The unit of a quantity-valued parameter read from metadata, whose value the transform holds as a plain number
# (degrees for an angle); `tag`, `unit_tag`, and `datatype` reproduce the stored quantity node. It applies only when
# the transform is written back with the tag family it was read from (`transform_tag`).
struct _ParameterUnit
    transform_tag::String
    property::String
    kind::Symbol
    unit::String
    tag::String
    unit_tag::String
    datatype::Union{Nothing,String}
end

# Properties of metadata nodes that GWCS.jl does not model, which writing restores: for the node at each path below the
# owner's node ("" for the owner's own), the tag family it was read with and its unmodeled properties.
const _PassThrough = Dict{String,Tuple{String,Dict{Any,Any}}}
# Shared by everything built without metadata, so such steps still compare equal; it is never mutated.
const _NO_PASS_THROUGH = _PassThrough()

"""
    AnnotatedTransform(transform; name=nothing, inputs=nothing, outputs=nothing)

A transform that evaluates `transform` unchanged and carries optional `name`, `inputs`, and `outputs` annotations, each
`nothing` when absent. Supplied names are nonempty, and the input and output names are each unique, with one per input
or output; otherwise construction throws. The inverse swaps the input and output names. [`compose`](@ref) and [`parallel`](@ref) never fuse an annotated transform into its neighbors.
[`annotate_transform`](@ref) changes some annotations and keeps the rest.

Python equivalent: the `name`, `inputs`, and `outputs` of `astropy.modeling.Model`.
"""
struct AnnotatedTransform{Nin,Nout,T<:AbstractTransform{Nin,Nout},N,IN,OUT} <: AbstractTransform{Nin,Nout}
    transform::T
    name::N
    inputs::IN
    outputs::OUT
    # Concrete whatever the units, so files that differ only in parameter units share one WCS type.
    parameter_units::Vector{_ParameterUnit}
    pass_through::_PassThrough
end

function _axis_names(names, count, boundary)
    names === nothing && return nothing
    values = Tuple(String.(names))
    length(values) == count || throw(DimensionMismatch("$boundary names have length $(length(values)), expected $count"))
    all(!isempty, values) || throw(ArgumentError("$boundary names must not be empty"))
    length(unique(values)) == count || throw(ArgumentError("$boundary names must be unique"))
    return values
end

_parameter_units(transform::AnnotatedTransform) = transform.parameter_units
_parameter_units(transform::AbstractTransform) = _ParameterUnit[]
_pass_through(transform::AnnotatedTransform) = transform.pass_through
_pass_through(transform::AbstractTransform) = _NO_PASS_THROUGH

function AnnotatedTransform(@nospecialize(transform::AbstractTransform{Nin,Nout}); name=nothing, inputs=nothing, outputs=nothing, parameter_units=_parameter_units(transform), pass_through=_pass_through(transform)) where {Nin,Nout}
    transform isa AnnotatedTransform && return annotate_transform(transform; name, inputs, outputs, parameter_units, pass_through)
    transform_name_value = name === nothing ? nothing : String(name)
    transform_name_value === nothing || !isempty(transform_name_value) || throw(ArgumentError("transform name must not be empty"))
    input_name_values = _axis_names(inputs, Nin, "input")
    output_name_values = _axis_names(outputs, Nout, "output")
    return AnnotatedTransform{Nin,Nout,typeof(transform),typeof(transform_name_value),typeof(input_name_values),typeof(output_name_values)}(transform, transform_name_value, input_name_values, output_name_values, parameter_units, pass_through)
end

"""
    annotate_transform(transform; name=transform_name(transform), inputs=input_names(transform),
                       outputs=output_names(transform)) -> AnnotatedTransform

Return `transform` as an [`AnnotatedTransform`](@ref) with the annotations given. An omitted keyword keeps the existing
[`transform_name`](@ref), [`input_names`](@ref), or [`output_names`](@ref), and `nothing` removes it.
"""
function annotate_transform(@nospecialize(transform::AbstractTransform); name=transform_name(transform), inputs=input_names(transform), outputs=output_names(transform), parameter_units=_parameter_units(transform), pass_through=_pass_through(transform))
    parent = transform isa AnnotatedTransform ? transform.transform : transform
    return AnnotatedTransform(parent; name, inputs, outputs, parameter_units, pass_through)
end

@inline (transform::AnnotatedTransform)(coords::StaticVector) = transform.transform(coords)

transform_name(transform::AnnotatedTransform) = transform.name
input_names(transform::AnnotatedTransform) = transform.inputs
output_names(transform::AnnotatedTransform) = transform.outputs

function inverse_transform(transform::AnnotatedTransform)
    inverse = _inverse_transform_checked(transform.transform)
    inverse isa NoInverse && return inverse
    return AnnotatedTransform(inverse; name=transform.name, inputs=transform.outputs, outputs=transform.inputs)
end

function compile_transform(transform::AnnotatedTransform)
    compiled = _compile_transform_checked(transform.transform)
    return AnnotatedTransform(compiled; name=transform.name, inputs=transform.inputs, outputs=transform.outputs, parameter_units=transform.parameter_units, pass_through=transform.pass_through)
end

_numeric_transform(transform::AnnotatedTransform) = _numeric_transform(transform.transform)

transform_jacobian(transform::AnnotatedTransform, coords::StaticVector) = _transform_jacobian_checked(transform.transform, coords)
