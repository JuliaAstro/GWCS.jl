"""A flattened sequence of transforms evaluated from left to right."""
struct ComposedTransform{Nin,Nout,Ts<:Tuple} <: AbstractTransform{Nin,Nout}
    transforms::Ts
end

function ComposedTransform(transforms::Tuple)
    isempty(transforms) && throw(ArgumentError("a composition requires at least one transform"))
    all(transform -> transform isa AbstractTransform, transforms) || throw(ArgumentError("all composed values must be transforms"))
    for i in 1:(length(transforms) - 1)
        output_dimension(transforms[i]) == input_dimension(transforms[i + 1]) || throw(DimensionMismatch("transform $i has $(output_dimension(transforms[i])) outputs but transform $(i + 1) has $(input_dimension(transforms[i + 1])) inputs"))
    end
    Nin = input_dimension(first(transforms))
    Nout = output_dimension(last(transforms))
    return ComposedTransform{Nin,Nout,typeof(transforms)}(transforms)
end

@inline _apply_composed(::Tuple{}, coords) = coords

@inline function _apply_composed(transforms::Tuple, coords)
    return _apply_composed(Base.tail(transforms), first(transforms)(coords))
end

@inline function (transform::ComposedTransform)(coords::StaticVector)
    return _apply_composed(transform.transforms, coords)
end

function _append_composed!(result, transform::ComposedTransform)
    foreach(item -> _append_composed!(result, item), transform.transforms)
    return result
end

function _append_composed!(result, transform::AbstractTransform)
    push!(result, transform)
    return result
end

_fuse_pair(first::Shift, second::Shift) = Shift(first.offset + second.offset)
_fuse_pair(first::Scale, second::Scale) = Scale(first.factor * second.factor)
_fuse_pair(first::Shift, second::Scale) = AffineTransform(SMatrix{1,1}(second.factor), SVector(second.factor * first.offset))
_fuse_pair(first::Scale, second::Shift) = AffineTransform(SMatrix{1,1}(first.factor), SVector(second.offset))

function _fuse_pair(first::AffineTransform{Nin,Nmid}, second::AffineTransform{Nmid,Nout}) where {Nin,Nmid,Nout}
    return AffineTransform(second.matrix * first.matrix, second.matrix * first.offset + second.offset)
end

_fuse_pair(first, second) = nothing

function compose(transforms::AbstractTransform...)
    isempty(transforms) && throw(ArgumentError("a composition requires at least one transform"))
    flattened = Any[]
    foreach(transform -> _append_composed!(flattened, transform), transforms)
    for i in 1:(length(flattened) - 1)
        output_dimension(flattened[i]) == input_dimension(flattened[i + 1]) || throw(DimensionMismatch("transform $i has $(output_dimension(flattened[i])) outputs but transform $(i + 1) has $(input_dimension(flattened[i + 1])) inputs"))
    end
    filter!(transform -> !(transform isa IdentityTransform), flattened)
    if isempty(flattened)
        N = input_dimension(first(transforms))
        return IdentityTransform{N}()
    end

    fused = Any[]
    for transform in flattened
        if !isempty(fused)
            replacement = _fuse_pair(last(fused), transform)
            if replacement !== nothing
                fused[end] = replacement
                continue
            end
        end
        push!(fused, transform)
    end

    length(fused) == 1 && return only(fused)
    return ComposedTransform(Tuple(fused))
end

Base.:|(first::AbstractTransform, second::AbstractTransform) = compose(first, second)

function inverse_transform(transform::ComposedTransform)
    inverses = map(_inverse_transform_checked, reverse(transform.transforms))
    any(inverse -> inverse isa NoInverse, inverses) && return NoInverse()
    return compose(inverses...)
end

@inline _inverse_initialize_composed(::Tuple{}, coords, used) = (coords, used)

@inline function _inverse_initialize_composed(transforms::Tuple, coords, used)
    initialized, local_used = _inverse_initialize(last(transforms), coords)
    initialized === nothing && return (nothing, false)
    return _inverse_initialize_composed(Base.front(transforms), initialized, used || local_used)
end

@inline _inverse_initialize(transform::ComposedTransform{N,N}, coords::StaticVector{N}) where {N} = _inverse_initialize_composed(transform.transforms, coords, false)

function compile_transform(transform::ComposedTransform)
    return compose(map(_compile_transform_checked, transform.transforms)...)
end

_numeric_transform(transform::ComposedTransform) = compose(map(_numeric_transform, transform.transforms)...)

input_names(transform::ComposedTransform) = input_names(first(transform.transforms))
output_names(transform::ComposedTransform) = output_names(last(transform.transforms))

@inline _apply_composed_jacobian(::Tuple{}, coords, jacobian) = (coords, jacobian)

@inline function _apply_composed_jacobian(transforms::Tuple, coords, jacobian)
    transform = first(transforms)
    local_jacobian = _transform_jacobian_checked(transform, coords)
    outputs = transform(coords)
    return _apply_composed_jacobian(Base.tail(transforms), outputs, local_jacobian * jacobian)
end

function transform_jacobian(transform::ComposedTransform{Nin}, coords::StaticVector{Nin}) where {Nin}
    T = _evaluation_type(eltype(coords))
    inputs = SVector{Nin,T}(coords)
    _, jacobian = _apply_composed_jacobian(transform.transforms, inputs, SMatrix{Nin,Nin,T}(I))
    return jacobian
end

"""Independent transforms evaluated over consecutive input slices."""
struct ParallelTransform{Nin,Nout,Ts<:Tuple} <: AbstractTransform{Nin,Nout}
    transforms::Ts
end

function ParallelTransform(transforms::Tuple)
    isempty(transforms) && throw(ArgumentError("a parallel transform requires at least one transform"))
    all(transform -> transform isa AbstractTransform, transforms) || throw(ArgumentError("all parallel values must be transforms"))
    Nin = sum(input_dimension, transforms)
    Nout = sum(output_dimension, transforms)
    return ParallelTransform{Nin,Nout,typeof(transforms)}(transforms)
end

@inline _apply_parallel(::Tuple{}, coords, input_offset) = ()

@inline function _apply_parallel(transforms::Tuple, coords, input_offset)
    transform = first(transforms)
    Nin = input_dimension(transform)
    inputs = SVector{Nin}(ntuple(i -> coords[input_offset + i], Nin))
    outputs = Tuple(transform(inputs))
    return (outputs..., _apply_parallel(Base.tail(transforms), coords, input_offset + Nin)...)
end

@inline function (transform::ParallelTransform{Nin,Nout})(coords::StaticVector{Nin}) where {Nin,Nout}
    return SVector{Nout}(_apply_parallel(transform.transforms, coords, 0))
end

function _append_parallel!(result, transform::ParallelTransform)
    foreach(item -> _append_parallel!(result, item), transform.transforms)
    return result
end

function _append_parallel!(result, transform::AbstractTransform)
    push!(result, transform)
    return result
end

_affine_form(transform::Shift) = (SMatrix{1,1}(one(transform.offset)), SVector(transform.offset))
_affine_form(transform::Scale) = (SMatrix{1,1}(transform.factor), SVector(zero(transform.factor)))
_affine_form(transform::AffineTransform) = (transform.matrix, transform.offset)
_affine_form(transform) = nothing

function _fuse_parallel_affine(transforms)
    forms = map(_affine_form, transforms)
    any(isnothing, forms) && return nothing
    Nin = sum(form -> size(form[1], 2), forms)
    Nout = sum(form -> size(form[1], 1), forms)
    T = promote_type(map(form -> eltype(form[1]), forms)...)
    matrix = MMatrix{Nout,Nin,T}(zeros(T, Nout, Nin))
    offset = MVector{Nout,T}(undef)
    input_offset = 0
    output_offset = 0
    for (block, block_offset) in forms
        nout, nin = size(block)
        matrix[(output_offset + 1):(output_offset + nout), (input_offset + 1):(input_offset + nin)] .= block
        offset[(output_offset + 1):(output_offset + nout)] .= block_offset
        input_offset += nin
        output_offset += nout
    end
    return AffineTransform(SMatrix(matrix), SVector(offset))
end

function parallel(transforms::AbstractTransform...)
    flattened = Any[]
    foreach(transform -> _append_parallel!(flattened, transform), transforms)
    isempty(flattened) && throw(ArgumentError("a parallel transform requires at least one transform"))
    fused = _fuse_parallel_affine(Tuple(flattened))
    fused !== nothing && return fused
    length(flattened) == 1 && return only(flattened)
    return ParallelTransform(Tuple(flattened))
end

Base.:&(first::AbstractTransform, second::AbstractTransform) = parallel(first, second)

function inverse_transform(transform::ParallelTransform)
    inverses = map(_inverse_transform_checked, transform.transforms)
    any(inverse -> inverse isa NoInverse, inverses) && return NoInverse()
    return parallel(inverses...)
end

@inline _inverse_initialize_parallel(::Tuple{}, coords, input_offset, output_offset) = ((), false)

@inline function _inverse_initialize_parallel(transforms::Tuple, coords, input_offset, output_offset)
    transform = first(transforms)
    Nout = output_dimension(transform)
    inputs = SVector{Nout}(ntuple(i -> coords[output_offset + i], Nout))
    initialized, used = _inverse_initialize(transform, inputs)
    initialized === nothing && return (nothing, false)
    tail, tail_used = _inverse_initialize_parallel(Base.tail(transforms), coords, input_offset + input_dimension(transform), output_offset + Nout)
    tail === nothing && return (nothing, false)
    return ((Tuple(initialized)..., tail...), used || tail_used)
end

@inline function _inverse_initialize(transform::ParallelTransform{Nin,Nout}, coords::StaticVector{Nout}) where {Nin,Nout}
    outputs, used = _inverse_initialize_parallel(transform.transforms, coords, 0, 0)
    return (outputs === nothing ? nothing : SVector{Nin}(outputs), used)
end

function compile_transform(transform::ParallelTransform)
    return parallel(map(_compile_transform_checked, transform.transforms)...)
end

_numeric_transform(transform::ParallelTransform) = parallel(map(_numeric_transform, transform.transforms)...)

function _parallel_names(transforms, accessor)
    names = map(accessor, transforms)
    all(name -> name !== nothing, names) || return nothing
    return Tuple(Iterators.flatten(names))
end

input_names(transform::ParallelTransform) = _parallel_names(transform.transforms, input_names)
output_names(transform::ParallelTransform) = _parallel_names(transform.transforms, output_names)

@inline _parallel_jacobians(::Tuple{}, coords, input_offset) = ()

@inline function _parallel_jacobians(transforms::Tuple, coords, input_offset)
    transform = first(transforms)
    Nin = input_dimension(transform)
    inputs = SVector{Nin}(ntuple(i -> coords[input_offset + i], Nin))
    jacobian = _transform_jacobian_checked(transform, inputs)
    return (jacobian, _parallel_jacobians(Base.tail(transforms), coords, input_offset + Nin)...)
end

@inline _parallel_jacobian_entry(::Tuple{}, row, column, input_offset, output_offset, ::Type{T}) where {T} = zero(T)

@inline function _parallel_jacobian_entry(jacobians::Tuple, row, column, input_offset, output_offset, ::Type{T}) where {T}
    jacobian = first(jacobians)
    Nout, Nin = size(jacobian)
    if output_offset < row <= output_offset + Nout && input_offset < column <= input_offset + Nin
        return T(jacobian[row - output_offset, column - input_offset])
    end
    return _parallel_jacobian_entry(Base.tail(jacobians), row, column, input_offset + Nin, output_offset + Nout, T)
end

function transform_jacobian(transform::ParallelTransform{Nin,Nout}, coords::StaticVector{Nin}) where {Nin,Nout}
    Tcoords = _evaluation_type(eltype(coords))
    inputs = SVector{Nin,Tcoords}(coords)
    jacobians = _parallel_jacobians(transform.transforms, inputs, 0)
    T = promote_type(map(eltype, jacobians)...)
    values = ntuple(Val(Nin * Nout)) do index
        row = mod1(index, Nout)
        column = div(index - 1, Nout) + 1
        _parallel_jacobian_entry(jacobians, row, column, 0, 0, T)
    end
    return SMatrix{Nout,Nin,T}(values)
end

abstract type _ArithmeticOperator end

struct _AddOperator <: _ArithmeticOperator end
struct _SubtractOperator <: _ArithmeticOperator end
struct _MultiplyOperator <: _ArithmeticOperator end
struct _DivideOperator <: _ArithmeticOperator end
struct _PowerOperator <: _ArithmeticOperator end

function _arithmetic_operator(operator)
    operator === (+) && return _AddOperator()
    operator === (-) && return _SubtractOperator()
    operator === (*) && return _MultiplyOperator()
    operator === (/) && return _DivideOperator()
    operator === (^) && return _PowerOperator()
    throw(ArgumentError("unsupported elementwise transform operator $operator"))
end

_apply_arithmetic(::_AddOperator, left, right) = left + right
_apply_arithmetic(::_SubtractOperator, left, right) = left - right
_apply_arithmetic(::_MultiplyOperator, left, right) = left * right
_apply_arithmetic(::_DivideOperator, left, right) = left / right
_apply_arithmetic(::_PowerOperator, left, right) = left ^ right

"""Two transforms evaluated on the same inputs and combined elementwise."""
struct ArithmeticTransform{Nin,Nout,O<:_ArithmeticOperator,L<:AbstractTransform{Nin,Nout},R<:AbstractTransform{Nin,Nout}} <: AbstractTransform{Nin,Nout}
    operator::O
    left::L
    right::R
end

function _arithmetic_transform(operator::_ArithmeticOperator, left::AbstractTransform, right::AbstractTransform)
    input_dimension(left) == input_dimension(right) || throw(DimensionMismatch("arithmetic operands have input dimensions $(input_dimension(left)) and $(input_dimension(right))"))
    output_dimension(left) == output_dimension(right) || throw(DimensionMismatch("arithmetic operands have output dimensions $(output_dimension(left)) and $(output_dimension(right))"))
    Nin = input_dimension(left)
    Nout = output_dimension(left)
    return ArithmeticTransform{Nin,Nout,typeof(operator),typeof(left),typeof(right)}(operator, left, right)
end

ArithmeticTransform(operator, left::AbstractTransform, right::AbstractTransform) = _arithmetic_transform(_arithmetic_operator(operator), left, right)

"""
    elementwise(operator, left, right)

Evaluate `left` and `right` on the same inputs and combine corresponding
outputs with `operator`.
"""
elementwise(operator, left::AbstractTransform, right::AbstractTransform) = ArithmeticTransform(operator, left, right)

Base.:+(left::AbstractTransform, right::AbstractTransform) = elementwise(+, left, right)
Base.:-(left::AbstractTransform, right::AbstractTransform) = elementwise(-, left, right)
Base.:*(left::AbstractTransform, right::AbstractTransform) = elementwise(*, left, right)
Base.:/(left::AbstractTransform, right::AbstractTransform) = elementwise(/, left, right)
Base.:^(left::AbstractTransform, right::AbstractTransform) = elementwise(^, left, right)

@inline function (transform::ArithmeticTransform{Nin,Nout})(coords::StaticVector{Nin}) where {Nin,Nout}
    left = transform.left(coords)
    right = transform.right(coords)
    return SVector{Nout}(ntuple(i -> _apply_arithmetic(transform.operator, left[i], right[i]), Nout))
end

function compile_transform(transform::ArithmeticTransform)
    left = _compile_transform_checked(transform.left)
    right = _compile_transform_checked(transform.right)
    return _arithmetic_transform(transform.operator, left, right)
end

_numeric_transform(transform::ArithmeticTransform) = _arithmetic_transform(transform.operator, _numeric_transform(transform.left), _numeric_transform(transform.right))

function _arithmetic_names(left, right, accessor)
    left_names = accessor(left)
    right_names = accessor(right)
    left_names === nothing && return right_names
    right_names === nothing && return left_names
    return left_names == right_names ? left_names : nothing
end

input_names(transform::ArithmeticTransform) = _arithmetic_names(transform.left, transform.right, input_names)
output_names(transform::ArithmeticTransform) = _arithmetic_names(transform.left, transform.right, output_names)

_arithmetic_jacobian_entry(::_AddOperator, left, right, left_derivative, right_derivative) = left_derivative + right_derivative
_arithmetic_jacobian_entry(::_SubtractOperator, left, right, left_derivative, right_derivative) = left_derivative - right_derivative
_arithmetic_jacobian_entry(::_MultiplyOperator, left, right, left_derivative, right_derivative) = right * left_derivative + left * right_derivative
_arithmetic_jacobian_entry(::_DivideOperator, left, right, left_derivative, right_derivative) = (right * left_derivative - left * right_derivative) / right^2

function _arithmetic_jacobian_entry(::_PowerOperator, left, right, left_derivative, right_derivative)
    value = left ^ right
    left_term = iszero(left_derivative) ? zero(value * left_derivative) : right * left ^ (right - one(right)) * left_derivative
    right_term = iszero(right_derivative) ? zero(value * right_derivative) : value * log(left) * right_derivative
    return left_term + right_term
end

function transform_jacobian(transform::ArithmeticTransform{Nin,Nout}, coords::StaticVector{Nin}) where {Nin,Nout}
    T = _evaluation_type(eltype(coords))
    inputs = SVector{Nin,T}(coords)
    left = transform.left(inputs)
    right = transform.right(inputs)
    left_jacobian = _transform_jacobian_checked(transform.left, inputs)
    right_jacobian = _transform_jacobian_checked(transform.right, inputs)
    values = ntuple(Val(Nin * Nout)) do index
        row = mod1(index, Nout)
        _arithmetic_jacobian_entry(transform.operator, left[row], right[row], left_jacobian[index], right_jacobian[index])
    end
    return SMatrix{Nout,Nin}(values)
end
