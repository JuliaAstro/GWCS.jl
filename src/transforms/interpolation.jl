"""A construction-time-copied multilinear lookup table over selected input axes."""
struct FITSLookupTable{N,Nin,T<:AbstractFloat,A<:AbstractArray{T,N},X<:NTuple{N,Int}} <: AbstractTransform{Nin,1}
    data::A
    axes::X
    crpix::SVector{N,T}
    crval::SVector{N,T}
    cdelt::SVector{N,T}
end

function FITSLookupTable(data::AbstractArray{<:Any,N}, axes::NTuple{N,Int}, ::Val{Nin}=Val(maximum(axes)); crpix=ntuple(_ -> 1, N), crval=ntuple(_ -> 0, N), cdelt=ntuple(_ -> 1, N)) where {N,Nin}
    N >= 1 || throw(ArgumentError("lookup table must have at least one dimension"))
    isempty(data) && throw(ArgumentError("lookup table data must not be empty"))
    allunique(axes) || throw(ArgumentError("lookup table axes must be distinct"))
    all(>(0), axes) || throw(ArgumentError("lookup table axes must be positive"))
    values = float.(data)
    T = eltype(values)
    reference_pixel = SVector{N,T}(crpix)
    reference_value = SVector{N,T}(crval)
    increment = SVector{N,T}(cdelt)
    all(isfinite, reference_pixel) && all(isfinite, reference_value) && all(isfinite, increment) || throw(ArgumentError("lookup table coordinates must be finite"))
    all(value -> !iszero(value), increment) || throw(ArgumentError("lookup table CDELT values must be nonzero"))
    all(<=(Nin), axes) || throw(ArgumentError("lookup table axis exceeds input dimension $Nin"))
    return FITSLookupTable{N,Nin,T,typeof(values),typeof(axes)}(values, axes, reference_pixel, reference_value, increment)
end

BilinearLookupTable(data::AbstractMatrix; kwargs...) = FITSLookupTable(data, (1, 2); kwargs...)

@inline function _lookup_cell(index, length)
    length == 1 && return 1, 1, zero(index)
    coordinate = clamp(index, one(index), typeof(index)(length))
    lower = min(floor(Int, coordinate), length - 1)
    return lower, lower + 1, coordinate - lower
end

@inline function (table::FITSLookupTable{1,Nin})(coords::StaticVector{Nin}) where {Nin}
    T = promote_type(eltype(coords), eltype(table.data))
    coordinate = (T(coords[table.axes[1]]) - T(table.crval[1])) / T(table.cdelt[1]) + T(table.crpix[1])
    isnan(coordinate) && return SVector(T(NaN))
    lower, upper, fraction = _lookup_cell(coordinate, length(table.data))
    return SVector(muladd(fraction, T(table.data[upper]) - T(table.data[lower]), T(table.data[lower])))
end

@inline function (table::FITSLookupTable{2,Nin})(coords::StaticVector{Nin}) where {Nin}
    T = promote_type(eltype(coords), eltype(table.data))
    x = (T(coords[table.axes[1]]) - T(table.crval[1])) / T(table.cdelt[1]) + T(table.crpix[1])
    y = (T(coords[table.axes[2]]) - T(table.crval[2])) / T(table.cdelt[2]) + T(table.crpix[2])
    (isnan(x) || isnan(y)) && return SVector(T(NaN))
    x0, x1, tx = _lookup_cell(x, size(table.data, 1))
    y0, y1, ty = _lookup_cell(y, size(table.data, 2))
    v00, v10 = T(table.data[x0, y0]), T(table.data[x1, y0])
    v01, v11 = T(table.data[x0, y1]), T(table.data[x1, y1])
    value = (1 - tx) * (1 - ty) * v00 + tx * (1 - ty) * v10 + (1 - tx) * ty * v01 + tx * ty * v11
    return SVector(value)
end

struct FITSLookupDistortion{N,D,S,C} <: AbstractTransform{N,N}
    detector_to_image::D
    sip::S
    cpdis::C
end

function FITSLookupDistortion(::Val{N}, detector_to_image, sip, cpdis) where {N}
    N >= 2 || throw(ArgumentError("FITS lookup distortions require at least two pixel axes"))
    length(detector_to_image) == 2 && length(cpdis) == 2 || throw(ArgumentError("FITS lookup distortion stages require two axis tables"))
    input_dimension(sip) == N && output_dimension(sip) == N || throw(DimensionMismatch("SIP transform must map $N axes to $N axes"))
    return FITSLookupDistortion{N,typeof(detector_to_image),typeof(sip),typeof(cpdis)}(detector_to_image, sip, cpdis)
end

@inline function _lookup_offsets(tables, coords::StaticVector{N,T}) where {N,T}
    x = tables[1] === nothing ? zero(T) : T(only(tables[1](coords)))
    y = tables[2] === nothing ? zero(T) : T(only(tables[2](coords)))
    return SVector{N,T}(ntuple(i -> i == 1 ? x : i == 2 ? y : zero(T), Val(N)))
end

@inline function (transform::FITSLookupDistortion{N})(coords::StaticVector{N}) where {N}
    detector = coords + _lookup_offsets(transform.detector_to_image, coords)
    corrected = transform.sip(detector)
    return corrected + _lookup_offsets(transform.cpdis, detector)
end

struct _FITSPriorDistortion{N,D,S,C,P} <: AbstractTransform{N,N}
    detector_to_image::D
    sip::S
    cpdis::C
    tpd::P
end

function _FITSPriorDistortion(::Val{N}, detector_to_image, sip, cpdis, tpd::_FITSTPDTransform{N}) where {N}
    length(detector_to_image) == 2 && length(cpdis) == 2 || throw(DimensionMismatch("FITS lookup distortion stages require two axis tables"))
    has_lookup = any(table -> table !== nothing, (detector_to_image..., cpdis...))
    has_lookup && N < 2 && throw(ArgumentError("FITS lookup distortions require at least two pixel axes"))
    input_dimension(sip) == N && output_dimension(sip) == N || throw(DimensionMismatch("SIP transform must map $N axes to $N axes"))
    return _FITSPriorDistortion{N,typeof(detector_to_image),typeof(sip),typeof(cpdis),typeof(tpd)}(detector_to_image, sip, cpdis, tpd)
end

@inline _fits_prior_outputs(::Tuple{}, ::Tuple{}, ::Tuple{}) = ()

@inline function _fits_prior_outputs(models::Tuple, lookup::Tuple, tpd::Tuple)
    output = first(models) === nothing ? first(lookup) : first(tpd)
    return (output, _fits_prior_outputs(Base.tail(models), Base.tail(lookup), Base.tail(tpd))...)
end

@inline function (transform::_FITSPriorDistortion{N})(coords::StaticVector{N}) where {N}
    detector = coords + _lookup_offsets(transform.detector_to_image, coords)
    corrected = transform.sip(detector)
    lookup_corrected = corrected + _lookup_offsets(transform.cpdis, detector)
    tpd_corrected = transform.tpd(corrected)
    outputs = _fits_prior_outputs(transform.tpd.models, Tuple(lookup_corrected), Tuple(tpd_corrected))
    return SVector{N}(outputs)
end

@inline function _inverse_initialize(transform::_FITSPriorDistortion{N}, coords::StaticVector{N}) where {N}
    tpd_initialized, tpd_used = _inverse_initialize(transform.tpd, coords)
    sip_initialized, sip_used = _inverse_initialize(transform.sip, tpd_initialized)
    return (sip_initialized, tpd_used || sip_used)
end

@inline function _inverse_initialize(transform::FITSLookupDistortion{N}, coords::StaticVector{N}) where {N}
    return _inverse_initialize(transform.sip, coords)
end

struct FITSTabularTable{M,A,I}
    axes::NTuple{M,Int}
    coordinates::A
    indices::I

    function FITSTabularTable{M,A,I}(axes, coordinates, indices) where {M,A,I}
        return new{M,A,I}(axes, coordinates, indices)
    end
end

function FITSTabularTable(axes::NTuple{M,Int}, coordinates, indices) where {M}
    M >= 1 || throw(ArgumentError("TAB tables must control at least one axis"))
    allunique(axes) || throw(ArgumentError("TAB table axes must be distinct"))
    values = Array{Float64}(coordinates)
    ndims(values) == M + 1 || throw(ArgumentError("TAB coordinate array has $(ndims(values) - 1) table axes, expected $M"))
    size(values, 1) == M || throw(ArgumentError("TAB coordinate array must have a leading component dimension of length $M"))
    all(isfinite, values) || throw(ArgumentError("TAB coordinate array values must be finite"))
    index_values = ntuple(M) do axis
        vector = Float64.(vec(indices[axis]))
        length(vector) == size(values, axis + 1) || throw(ArgumentError("TAB index vector $axis has length $(length(vector)), expected $(size(values, axis + 1))"))
        all(isfinite, vector) || throw(ArgumentError("TAB index vector $axis values must be finite"))
        if length(vector) > 1
            increasing = all(i -> vector[i] <= vector[i + 1], 1:(length(vector) - 1))
            decreasing = all(i -> vector[i] >= vector[i + 1], 1:(length(vector) - 1))
            (increasing || decreasing) && vector[1] != vector[end] || throw(ArgumentError("TAB index vector $axis must be monotonic"))
        end
        vector
    end
    return FITSTabularTable{M,typeof(values),typeof(index_values)}(axes, values, index_values)
end

# First interval `index[i:i+1]` containing `value` with distinct endpoints; runs of equal
# index values are skipped because they define no interpolation interval.
# Takes an ordering object because `searchsortedfirst(...; rev)` with a runtime Bool allocates on Julia 1.10.
@inline function _tabular_interval(index, value, order)
    upper = searchsortedfirst(index, value, order)
    return upper > 1 ? upper - 1 : searchsortedlast(index, value, order)
end

function _tabular_fractional_index(index, value)
    T = promote_type(eltype(index), typeof(float(value)))
    if length(index) == 1
        abs(value - only(index)) <= T(0.5) || throw(DomainError(value, "TAB coordinate is outside index-vector bounds"))
        return one(T)
    end
    increasing = index[end] > index[1]
    if increasing
        index[1] - (index[2] - index[1]) / 2 <= value <= index[end] + (index[end] - index[end - 1]) / 2 || throw(DomainError(value, "TAB coordinate is outside index-vector bounds"))
        position = value < index[1] ? 1 : value > index[end] ? length(index) - 1 : _tabular_interval(index, value, Base.Order.Forward)
    else
        index[end] - (index[end - 1] - index[end]) / 2 <= value <= index[1] + (index[1] - index[2]) / 2 || throw(DomainError(value, "TAB coordinate is outside index-vector bounds"))
        position = value > index[1] ? 1 : value < index[end] ? length(index) - 1 : _tabular_interval(index, value, Base.Order.Reverse)
    end
    1 <= position < length(index) || throw(DomainError(value, "TAB coordinate does not lie in a nondegenerate index-vector interval"))
    first, second = T(index[position]), T(index[position + 1])
    return T(position) + (value - first) / (second - first)
end

@inline function _tabular_cell(coordinate, length)
    length == 1 && return 1, 1, zero(coordinate)
    T = typeof(coordinate)
    T(0.5) <= coordinate <= T(length) + T(0.5) || throw(DomainError(coordinate, "TAB fractional index is outside interpolation bounds"))
    lower = floor(Int, coordinate)
    fraction = coordinate - lower
    lower == 0 && return 1, 2, fraction - 1
    lower == length && return length - 1, length, fraction + 1
    return lower, lower + 1, fraction
end

@inline function _tabular_index_value(index, coordinate)
    length(index) == 1 && return oftype(coordinate, only(index))
    lower, upper, weight = _tabular_cell(coordinate, length(index))
    return (1 - weight) * index[lower] + weight * index[upper]
end

@inline function _tabular_values(table::FITSTabularTable{M}, coords) where {M}
    any(axis -> !isfinite(coords[axis]), table.axes) && return ntuple(_ -> promote_type(eltype(coords), eltype(table.coordinates))(NaN), M)
    cells = ntuple(M) do axis
        fractional = _tabular_fractional_index(table.indices[axis], coords[table.axes[axis]])
        _tabular_cell(fractional, size(table.coordinates, axis + 1))
    end
    return ntuple(M) do component
        total = zero(coords[1])
        for corner in 0:(2^M - 1)
            indices = ntuple(M) do axis
                lower, upper, fraction = cells[axis]
                high = !iszero((corner >> (axis - 1)) & 1)
                high ? upper : lower
            end
            weight = prod(ntuple(M) do axis
                fraction = cells[axis][3]
                iszero((corner >> (axis - 1)) & 1) ? 1 - fraction : fraction
            end)
            total += weight * table.coordinates[component, indices...]
        end
        total
    end
end

struct FITSTabularTransform{N,T,S} <: AbstractTransform{N,N}
    tables::T
    scales::S
end

function FITSTabularTransform(::Val{N}, tables, scales) where {N}
    length(scales) == N || throw(DimensionMismatch("TAB output scaling has length $(length(scales)), expected $N"))
    used = Int[]
    for table in tables, axis in table.axes
        1 <= axis <= N || throw(ArgumentError("TAB table references axis $axis outside 1:$N"))
        axis in used && throw(ArgumentError("multiple TAB tables control axis $axis"))
        push!(used, axis)
    end
    return FITSTabularTransform{N,typeof(tables),typeof(Tuple(scales))}(tables, Tuple(scales))
end

# Recurse over the table tuple so tables with different dimensions stay type-stable.
@inline _tabular_apply(::Tuple{}, coords, result, scales) = result

@inline function _tabular_apply(tables::Tuple, coords, result, scales)
    table = first(tables)
    values = _tabular_values(table, coords)
    for component in eachindex(table.axes)
        axis = table.axes[component]
        result = Base.setindex(result, values[component] * scales[axis], axis)
    end
    return _tabular_apply(Base.tail(tables), coords, result, scales)
end

@inline (transform::FITSTabularTransform{N})(coords::StaticVector{N}) where {N} = _tabular_apply(transform.tables, coords, SVector(coords), transform.scales)

struct FITSTabularInverse{N,T,S} <: AbstractTransform{N,N}
    tables::T
    scales::S
end

function inverse_transform(transform::FITSTabularTransform{N}) where {N}
    all(table -> length(table.axes) == 1 && _tabular_monotonic(table.coordinates), transform.tables) || return NoInverse()
    return FITSTabularInverse{N,typeof(transform.tables),typeof(transform.scales)}(transform.tables, transform.scales)
end

_tabular_monotonic(values) = length(values) <= 1 || values[1] != values[end] &&
    (all(i -> values[i] <= values[i + 1], 1:(length(values) - 1)) || all(i -> values[i] >= values[i + 1], 1:(length(values) - 1)))

@inline function (transform::FITSTabularInverse{N})(coords::StaticVector{N}) where {N}
    result = SVector(coords)
    for table in transform.tables
        axis = only(table.axes)
        isfinite(coords[axis]) || (result = Base.setindex(result, NaN, axis); continue)
        fractional = _tabular_fractional_index(view(table.coordinates, 1, :), coords[axis] / transform.scales[axis])
        result = Base.setindex(result, _tabular_index_value(only(table.indices), fractional), axis)
    end
    return SVector(result)
end

inverse_transform(transform::FITSTabularInverse{N}) where {N} = FITSTabularTransform(Val(N), transform.tables, transform.scales)
