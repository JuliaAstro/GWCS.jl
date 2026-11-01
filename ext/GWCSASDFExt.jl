module GWCSASDFExt

using ASDF
using GWCS
using PrecompileTools: @compile_workload

GWCS._tag(node::ASDF.TaggedMapping) = node.tag
GWCS._tag(node::ASDF.TaggedSequence) = node.tag
GWCS._tag(node::ASDF.TaggedScalar) = node.tag
GWCS._unwrap(node::ASDF.TaggedScalar) = String(node)
GWCS._has_tag(::Union{ASDF.TaggedMapping,ASDF.TaggedSequence,ASDF.TaggedScalar}) = true

function GWCS._materialize_array(array::ASDF.NDArray)
    value = array[]
    return value isa AbstractMatrix ? permutedims(value) : value
end

# An unmodeled array keeps its element type and is written back as an ndarray.
GWCS._pass_through_value(array::Union{ASDF.NDArray,ASDF.ChunkedNDArray}) = GWCS._MetadataArray(_reverse_axes(array[]))
_reverse_axes(array) = permutedims(array, ndims(array):-1:1)

const _COMPRESSIONS = Dict(:zlib => ASDF.C_Zlib, :bzp2 => ASDF.C_Bzip2, :lz4 => ASDF.C_Lz4, :zstd => ASDF.C_Zstd)

function _array_policy(threshold, storage, compression)
    storage === :external && throw(ArgumentError("all_array_storage=:external is not supported because ASDF.jl cannot write external array blocks"))
    storage in (nothing, :internal, :inline) || throw(ArgumentError("all_array_storage must be nothing, :internal, or :inline, got $(repr(storage))"))
    compression === nothing || haskey(_COMPRESSIONS, compression) || throw(ArgumentError("all_array_compression must be nothing, :zlib, :bzp2, :lz4, or :zstd, got $(repr(compression))"))
    return (; threshold, storage, compression=compression === nothing ? nothing : _COMPRESSIONS[compression])
end

# Python asdf's precedence: `all_array_storage` replaces an array's own storage, and a threshold then chooses inline or
# block storage by size.
function _ndarray(array, policy; inline=false, compression=ASDF.C_None, lz4_layout=:block)
    policy.storage === nothing || (inline = policy.storage === :inline)
    policy.threshold === nothing || (inline = length(array) < policy.threshold)
    policy.compression === nothing || (compression = policy.compression)
    return ASDF.NDArrayWrapper(array; inline, compression, lz4_layout)
end

_storable(array) = haskey(ASDF.type_datatype_dict, eltype(array))

_asdf_metadata(value::GWCS.TaggedMapping, policy) = ASDF.TaggedMapping(value.tag, Dict{Any,Any}(_asdf_metadata(key, policy) => _asdf_metadata(item, policy) for (key, item) in value))
_asdf_metadata(value::GWCS.TaggedSequence, policy) = ASDF.TaggedSequence(value.tag, [_asdf_metadata(item, policy) for item in value.value])
_asdf_metadata(value::GWCS.TaggedScalar, policy) = ASDF.TaggedScalar(value.tag, String(value.value))
_asdf_metadata(value::AbstractDict, policy) = Dict{Any,Any}(_asdf_metadata(key, policy) => _asdf_metadata(item, policy) for (key, item) in value)
# ASDF.jl writes vectors as YAML lists and requires other arrays to be ndarrays.
_asdf_metadata(value::AbstractArray, policy) = ndims(value) >= 2 && _storable(value) ? _ndarray(value, policy) : map(item -> _asdf_metadata(item, policy), value)
_asdf_metadata(value::Tuple, policy) = map(item -> _asdf_metadata(item, policy), value)
# Labels ASDF.jl cannot store as an ndarray, such as strings, stay nested lists.
_asdf_metadata(value::GWCS._MetadataArray, policy) = _storable(value.array) ? _ndarray(_reverse_axes(value.array), policy) : _asdf_metadata(GWCS._nested_lists(value.array), policy)
_asdf_metadata(value::Union{ASDF.NDArray,ASDF.ChunkedNDArray}, policy) = _ndarray(value[], policy)
_asdf_metadata(value::ASDF.NDArrayWrapper, policy) = _ndarray(value.array, policy; value.inline, value.compression, value.lz4_layout)
_asdf_metadata(value::Union{GWCS.WCS,GWCS.Step,GWCS.AbstractTransform,GWCS.AbstractCoordinateFrame,GWCS.AbstractBoundingBox,GWCS.AbstractLabelMapper}, policy) = _asdf_metadata(GWCS._to_metadata(value, "metadata", GWCS.default_transform_registry()), policy)
_asdf_metadata(value, policy) = value

function GWCS.write_asdf(filename, wcs::GWCS.WCS; key="wcs", array_inline_threshold=nothing, all_array_storage=nothing, all_array_compression=nothing)
    return GWCS.write_asdf(filename, Dict(String(key) => wcs); array_inline_threshold, all_array_storage, all_array_compression)
end

function GWCS.write_asdf(filename, document::AbstractDict; array_inline_threshold::Union{Nothing,Integer}=nothing, all_array_storage=nothing, all_array_compression=nothing)
    return ASDF.write_file(filename, _asdf_metadata(document, _array_policy(array_inline_threshold, all_array_storage, all_array_compression)))
end

# Parses the main workload's pipeline as the OrderedDict-backed nodes `ASDF.load_file` builds, so the cached
# adapter and parser specializations match loaded files.
_ordered(value::ASDF.TaggedMapping) = ASDF.TaggedMapping(value.tag, ASDF.OrderedDict{Any,Any}(_ordered(k) => _ordered(v) for (k, v) in value.value))
_ordered(value::ASDF.TaggedSequence) = ASDF.TaggedSequence(value.tag, map(_ordered, value.value))
_ordered(value::AbstractDict) = ASDF.OrderedDict{Any,Any}(_ordered(k) => _ordered(v) for (k, v) in value)
_ordered(value::AbstractArray) = map(_ordered, value)
_ordered(value) = value
@compile_workload begin
    wcs = GWCS.parse_wcs(_ordered(_asdf_metadata(GWCS._load_tagged_yaml(GWCS._WORKLOAD_ASDF), _array_policy(nothing, nothing, nothing))))
    GWCS.world_to_pixel(wcs, GWCS.pixel_to_world(wcs, 100.0, 200.0))
end

end
