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

"""
    AnnotatedTransform(transform; name=nothing, inputs=nothing, outputs=nothing)
    annotate_transform(transform; name, inputs, outputs)

Attach an optional transform name and input/output axis names without changing
the wrapped transform's numerical behavior.
"""
struct AnnotatedTransform{Nin,Nout,T<:AbstractTransform{Nin,Nout},N,IN,OUT} <: AbstractTransform{Nin,Nout}
    transform::T
    name::N
    inputs::IN
    outputs::OUT
    # Concrete whatever the units, so files that differ only in parameter units share one WCS type.
    parameter_units::Vector{_ParameterUnit}
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

function AnnotatedTransform(@nospecialize(transform::AbstractTransform{Nin,Nout}); name=nothing, inputs=nothing, outputs=nothing, parameter_units=_parameter_units(transform)) where {Nin,Nout}
    transform isa AnnotatedTransform && return annotate_transform(transform; name, inputs, outputs, parameter_units)
    transform_name_value = name === nothing ? nothing : String(name)
    transform_name_value === nothing || !isempty(transform_name_value) || throw(ArgumentError("transform name must not be empty"))
    input_name_values = _axis_names(inputs, Nin, "input")
    output_name_values = _axis_names(outputs, Nout, "output")
    return AnnotatedTransform{Nin,Nout,typeof(transform),typeof(transform_name_value),typeof(input_name_values),typeof(output_name_values)}(transform, transform_name_value, input_name_values, output_name_values, parameter_units)
end

function annotate_transform(@nospecialize(transform::AbstractTransform); name=transform_name(transform), inputs=input_names(transform), outputs=output_names(transform), parameter_units=_parameter_units(transform))
    parent = transform isa AnnotatedTransform ? transform.transform : transform
    return AnnotatedTransform(parent; name, inputs, outputs, parameter_units)
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
    return AnnotatedTransform(compiled; name=transform.name, inputs=transform.inputs, outputs=transform.outputs, parameter_units=transform.parameter_units)
end

_numeric_transform(transform::AnnotatedTransform) = _numeric_transform(transform.transform)

transform_jacobian(transform::AnnotatedTransform, coords::StaticVector) = _transform_jacobian_checked(transform.transform, coords)
