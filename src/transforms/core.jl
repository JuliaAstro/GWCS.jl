"""
    AbstractTransform{Nin,Nout}

Supertype of every transform, with `Nin` inputs and `Nout` outputs fixed by the type. A transform `t` is callable:
`t(coords)` takes a `StaticVector` of length `Nin` and returns an `SVector` of length `Nout`, with `NaN` in every output
where the coordinate has no image.

A new transform subtypes `AbstractTransform`, defines the call, and extends the protocol where it can do better than the
defaults: [`inverse_transform`](@ref), [`transform_jacobian`](@ref), [`transform_dependencies`](@ref),
[`compile_transform`](@ref), and the names [`transform_name`](@ref), [`input_names`](@ref), and
[`output_names`](@ref). [`register_transform!`](@ref) and [`register_transform_serializer!`](@ref) let metadata read
and write it.

Python equivalent: `astropy.modeling.Model`.
"""
abstract type AbstractTransform{Nin,Nout} end

"""
    NoInverse()

The value [`inverse_transform`](@ref) and [`get_transform`](@ref) return when there is no inverse.
"""
struct NoInverse end

"""
    input_dimension(x) -> Int

Return the number of inputs of a transform, label mapper, or bounding box.

Python equivalent: `astropy.modeling.Model.n_inputs`.
"""
input_dimension(::AbstractTransform{Nin,Nout}) where {Nin,Nout} = Nin

"""
    output_dimension(transform) -> Int

Return the number of outputs of a transform.

Python equivalent: `astropy.modeling.Model.n_outputs`.
"""
output_dimension(::AbstractTransform{Nin,Nout}) where {Nin,Nout} = Nout

"""
    transform_name(transform) -> Union{String,Nothing}

Return the name of a transform, or `nothing` when it has none. [`AnnotatedTransform`](@ref) and metadata give names.

Python equivalent: `astropy.modeling.Model.name`.
"""
transform_name(::AbstractTransform) = nothing

"""
    input_names(transform) -> Union{NTuple{input_dimension,String},Nothing}

Return the names of a transform's inputs, or `nothing` when they are unnamed. [`AnnotatedTransform`](@ref) and metadata
give names, and [`fix_inputs`](@ref) accepts them.

Python equivalent: `astropy.modeling.Model.inputs`, which always has names.
"""
input_names(::AbstractTransform) = nothing

"""
    output_names(transform) -> Union{NTuple{output_dimension,String},Nothing}

Return the names of a transform's outputs, or `nothing` when they are unnamed.

Python equivalent: `astropy.modeling.Model.outputs`, which always has names.
"""
output_names(::AbstractTransform) = nothing

function _dense_dependencies(::Val{Nout}, ::Val{Nin}) where {Nin,Nout}
    return SMatrix{Nout,Nin,Bool}(ntuple(_ -> true, Nin * Nout))
end

"""
    transform_dependencies(transform) -> SMatrix{output_dimension,input_dimension,Bool}

Return which outputs of `transform` depend on which inputs: entry `(i, j)` is `true` when output `i` can change with
input `j`. A `false` entry holds for every input and parameter value, not only where a derivative happens to be zero.
The default assumes every output depends on every input; a transform with independent axes extends this method.

Python equivalent: `astropy.modeling.separable.separability_matrix`.
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
    inverse_transform(transform) -> Union{AbstractTransform,NoInverse}

Return the analytical or explicitly supplied inverse of `transform`, or [`NoInverse`](@ref)`()` when it has none. A
combination of transforms has an inverse only when each part does. A new transform extends this method to supply its
inverse.

Python equivalent: `astropy.modeling.Model.inverse`, which raises `NotImplementedError` when there is none.
"""
inverse_transform(::AbstractTransform) = NoInverse()

struct _ApproximateInverseTransform{N,F,R} <: AbstractTransform{N,N}
    forward::F
    reverse::R
end

_ApproximateInverseTransform(forward::AbstractTransform{N,N}, reverse::AbstractTransform{N,N}) where {N} = _ApproximateInverseTransform{N,typeof(forward),typeof(reverse)}(forward, reverse)
@inline (transform::_ApproximateInverseTransform)(coords::StaticVector) = transform.forward(coords)
@inline _inverse_initialize(transform::_ApproximateInverseTransform{N}, coords::StaticVector{N}) where {N} = (transform.reverse(coords), true)
transform_jacobian(transform::_ApproximateInverseTransform, coords::StaticVector) = _transform_jacobian_checked(transform.forward, coords)

function _inverse_initialize(transform::AbstractTransform{Nin,Nout}, coords::StaticVector{Nout}) where {Nin,Nout}
    inverse = _inverse_transform_checked(transform)
    inverse isa NoInverse || return (inverse(coords), false)
    return Nin == Nout ? (coords, false) : (nothing, false)
end

function _inverse_transform_checked(transform::AbstractTransform)
    inverse = inverse_transform(transform)
    inverse isa NoInverse && return inverse
    inverse isa AbstractTransform || throw(ArgumentError("inverse_transform(::$(typeof(transform))) returned $(typeof(inverse)), expected AbstractTransform or NoInverse"))
    input_dimension(inverse) == output_dimension(transform) && output_dimension(inverse) == input_dimension(transform) || throw(DimensionMismatch("inverse of $(typeof(transform)) has dimensions $(input_dimension(inverse)) -> $(output_dimension(inverse)), expected $(output_dimension(transform)) -> $(input_dimension(transform))"))
    return inverse
end

"""
    has_inverse(transform) -> Bool

Return whether [`inverse_transform`](@ref) gives an inverse for `transform` rather than [`NoInverse`](@ref).

Python equivalent: `astropy.modeling.Model.has_inverse`.
"""
has_inverse(transform::AbstractTransform) = !(_inverse_transform_checked(transform) isa NoInverse)

"""
    compile_transform(transform) -> AbstractTransform

Return the form of `transform` that a [`WCS`](@ref) evaluates, with the same dimensions. The default returns
`transform`; a transform built for convenience extends this method to return a faster equivalent.
"""
compile_transform(transform::AbstractTransform) = transform
_numeric_transform(transform::AbstractTransform) = transform

function _compile_transform_checked(@nospecialize(transform::AbstractTransform))
    compiled = compile_transform(transform)
    compiled isa AbstractTransform || throw(ArgumentError("compile_transform(::$(typeof(transform))) returned $(typeof(compiled)), expected AbstractTransform"))
    input_dimension(compiled) == input_dimension(transform) && output_dimension(compiled) == output_dimension(transform) || throw(DimensionMismatch("compiled form of $(typeof(transform)) has dimensions $(input_dimension(compiled)) -> $(output_dimension(compiled)), expected $(input_dimension(transform)) -> $(output_dimension(transform))"))
    return compiled
end

function _forwarddiff_jacobian(transform::AbstractTransform{Nin}, coords::StaticVector{Nin}) where {Nin}
    T = _evaluation_type(eltype(coords))
    return ForwardDiff.jacobian(transform, SVector{Nin,T}(coords))
end

"""
    transform_jacobian(transform, coords::StaticVector) -> SMatrix{output_dimension,input_dimension}

Return the Jacobian of `transform` at `coords`, the matrix of `∂output_i/∂input_j`. The default differentiates the
transform with ForwardDiff.jl. A transform extends this method to supply an analytical Jacobian, which every WCS
Jacobian then uses.
"""
transform_jacobian(transform::AbstractTransform, coords::StaticVector) = _forwarddiff_jacobian(transform, coords)

function _transform_jacobian_checked(transform::AbstractTransform{Nin,Nout}, coords::StaticVector{Nin}) where {Nin,Nout}
    jacobian = transform_jacobian(transform, coords)
    jacobian isa AbstractMatrix || throw(ArgumentError("transform_jacobian(::$(typeof(transform)), ::$(typeof(coords))) returned $(typeof(jacobian)), expected an AbstractMatrix"))
    size(jacobian) == (Nout, Nin) || throw(DimensionMismatch("Jacobian of $(typeof(transform)) has size $(size(jacobian)), expected ($Nout, $Nin)"))
    return SMatrix{Nout,Nin}(jacobian)
end

_value_and_jacobian(transform::AbstractTransform, coords::StaticVector) = (transform(coords), _transform_jacobian_checked(transform, coords))

_evaluation_type(::Type{T}) where {T<:Integer} = Float64
_evaluation_type(::Type{T}) where {T<:Rational} = Float64
_evaluation_type(::Type{T}) where {T} = T

"""
    WithInverse(forward, backward)

A transform that evaluates `forward` and uses `backward` as its explicit inverse. `backward` maps the outputs of
`forward` to its inputs. GWCS.jl does not check that it inverts `forward` exactly, so an approximate inverse such as a
fitted polynomial is allowed. Its inverse is `WithInverse(backward, forward)`.

Python equivalent: setting `astropy.modeling.Model.inverse`.
"""
struct WithInverse{Nin,Nout,F<:AbstractTransform{Nin,Nout},I<:AbstractTransform{Nout,Nin}} <: AbstractTransform{Nin,Nout}
    forward::F
    backward::I
end

@inline (transform::WithInverse)(coords::StaticVector) = transform.forward(coords)

inverse_transform(transform::WithInverse) = WithInverse(transform.backward, transform.forward)
compile_transform(transform::WithInverse) = WithInverse(_compile_transform_checked(transform.forward), _compile_transform_checked(transform.backward))
_numeric_transform(transform::WithInverse) = WithInverse(_numeric_transform(transform.forward), _numeric_transform(transform.backward))
transform_jacobian(transform::WithInverse, coords::StaticVector) = _transform_jacobian_checked(transform.forward, coords)
transform_name(transform::WithInverse) = transform_name(transform.forward)
input_names(transform::WithInverse) = input_names(transform.forward)
output_names(transform::WithInverse) = output_names(transform.forward)
