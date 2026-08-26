"""Identity transform in `N` dimensions."""
struct IdentityTransform{N} <: AbstractTransform{N,N} end

@inline (::IdentityTransform{N})(coords::StaticVector{N}) where {N} = coords
inverse_transform(transform::IdentityTransform) = transform
transform_jacobian(::IdentityTransform{N}, coords::StaticVector{N}) where {N} = SMatrix{N,N,_evaluation_type(eltype(coords))}(I)

"""Add a constant offset to one coordinate."""
struct Shift{T<:Real} <: AbstractTransform{1,1}
    offset::T
end

@inline function (transform::Shift)(coords::StaticVector{1})
    T = _evaluation_type(eltype(coords))
    x = T(coords[1])
    return SVector(x + T(transform.offset))
end

inverse_transform(transform::Shift) = Shift(-transform.offset)
transform_jacobian(::Shift, coords::StaticVector{1}) = SMatrix{1,1}(one(_evaluation_type(eltype(coords))))

"""Multiply one coordinate by a constant factor."""
struct Scale{T<:Real} <: AbstractTransform{1,1}
    factor::T
end

@inline function (transform::Scale)(coords::StaticVector{1})
    T = _evaluation_type(eltype(coords))
    x = T(coords[1])
    return SVector(x * T(transform.factor))
end

function inverse_transform(transform::Scale)
    iszero(transform.factor) && return NoInverse()
    return Scale(inv(transform.factor))
end

function transform_jacobian(transform::Scale, coords::StaticVector{1})
    T = _evaluation_type(eltype(coords))
    return SMatrix{1,1}(T(transform.factor))
end

"""Static affine transform `matrix * coords + offset`."""
struct AffineTransform{Nin,Nout,T,L} <: AbstractTransform{Nin,Nout}
    matrix::SMatrix{Nout,Nin,T,L}
    offset::SVector{Nout,T}
end

function AffineTransform(matrix::StaticMatrix{Nout,Nin}, offset::StaticVector{Nout}) where {Nin,Nout}
    T = promote_type(eltype(matrix), eltype(offset))
    return AffineTransform{Nin,Nout,T,Nin * Nout}(SMatrix{Nout,Nin,T}(matrix), SVector{Nout,T}(offset))
end

function AffineTransform(matrix::AbstractMatrix, offset::AbstractVector)
    size(matrix, 1) == length(offset) || throw(DimensionMismatch("matrix has $(size(matrix, 1)) outputs but offset has length $(length(offset))"))
    Nin = size(matrix, 2)
    Nout = size(matrix, 1)
    T = promote_type(eltype(matrix), eltype(offset))
    return AffineTransform(SMatrix{Nout,Nin,T}(matrix), SVector{Nout,T}(offset))
end

function AffineTransform(matrix::AbstractMatrix)
    return AffineTransform(matrix, zeros(eltype(matrix), size(matrix, 1)))
end

@inline function (transform::AffineTransform{Nin,Nout})(coords::StaticVector{Nin}) where {Nin,Nout}
    T = _evaluation_type(eltype(coords))
    return SMatrix{Nout,Nin,T}(transform.matrix) * SVector{Nin,T}(coords) + SVector{Nout,T}(transform.offset)
end

function inverse_transform(transform::AffineTransform{N,N}) where {N}
    iszero(det(transform.matrix)) && return NoInverse()
    matrix = inv(transform.matrix)
    return AffineTransform(matrix, -(matrix * transform.offset))
end

inverse_transform(::AffineTransform) = NoInverse()

function transform_jacobian(transform::AffineTransform{Nin,Nout}, coords::StaticVector{Nin}) where {Nin,Nout}
    T = _evaluation_type(eltype(coords))
    return SMatrix{Nout,Nin,T}(transform.matrix)
end

"""
    RemapAxes{Nin}(mapping)

Reorder, duplicate, or drop inputs. `mapping` uses Julia's one-based indices.
"""
struct RemapAxes{Nin,Nout,M<:NTuple{Nout,Int}} <: AbstractTransform{Nin,Nout}
    mapping::M
end

function RemapAxes{Nin}(mapping::NTuple{Nout,Int}) where {Nin,Nout}
    all(i -> 1 <= i <= Nin, mapping) || throw(ArgumentError("mapping indices must be between 1 and $Nin"))
    return RemapAxes{Nin,Nout,typeof(mapping)}(mapping)
end

RemapAxes{Nin}(mapping::AbstractVector{<:Integer}) where {Nin} = RemapAxes{Nin}(Tuple(Int.(mapping)))

@inline function (transform::RemapAxes{Nin,Nout})(coords::StaticVector{Nin}) where {Nin,Nout}
    return SVector{Nout}(ntuple(i -> coords[transform.mapping[i]], Nout))
end

function inverse_transform(transform::RemapAxes{Nin,Nout}) where {Nin,Nout}
    all(i -> i in transform.mapping, 1:Nin) || return NoInverse()
    mapping = ntuple(i -> findfirst(==(i), transform.mapping)::Int, Nin)
    return RemapAxes{Nout}(mapping)
end
