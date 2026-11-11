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

# `ASDF.write_file` and `save` write GWCS objects anywhere in a document, with their arrays in uncompressed blocks by default.
const _GWCSObject = Union{GWCS.WCS,GWCS.Step,GWCS.AbstractTransform,GWCS.AbstractCoordinateFrame,GWCS.AbstractBoundingBox,GWCS.AbstractLabelMapper}
ASDF.to_tree(value::_GWCSObject) = GWCS._to_metadata(value, "metadata", GWCS.default_transform_registry())
ASDF.to_tree(value::GWCS.TaggedMapping) = ASDF.TaggedMapping(value.tag, value.value)
ASDF.to_tree(value::GWCS.TaggedSequence) = ASDF.TaggedSequence(value.tag, value.value)
ASDF.to_tree(value::GWCS.TaggedScalar) = ASDF.TaggedScalar(value.tag, String(value.value))
ASDF.to_tree(value::GWCS._MetadataArray) = ASDF.NDArrayWrapper(_reverse_axes(value.array); compression=ASDF.C_None)

GWCS.write_asdf(filename, wcs::GWCS.WCS; key="wcs") = ASDF.write_file(filename, Dict(String(key) => wcs))

# Parses the main workload's pipeline as the OrderedDict-backed nodes `ASDF.load_file` builds, so the cached
# adapter and parser specializations match loaded files.
_ordered(value::GWCS.TaggedMapping) = ASDF.TaggedMapping(value.tag, ASDF.OrderedDict{Any,Any}(_ordered(k) => _ordered(v) for (k, v) in value.value))
_ordered(value::GWCS.TaggedSequence) = ASDF.TaggedSequence(value.tag, map(_ordered, value.value))
_ordered(value::GWCS.TaggedScalar) = ASDF.TaggedScalar(value.tag, String(value.value))
_ordered(value::AbstractDict) = ASDF.OrderedDict{Any,Any}(_ordered(k) => _ordered(v) for (k, v) in value)
_ordered(value::AbstractArray) = map(_ordered, value)
_ordered(value) = value
@compile_workload begin
    wcs = GWCS.parse_wcs(_ordered(GWCS._load_tagged_yaml(GWCS._WORKLOAD_ASDF)))
    GWCS.world_to_pixel(wcs, GWCS.pixel_to_world(wcs, 100.0, 200.0))
end

end
