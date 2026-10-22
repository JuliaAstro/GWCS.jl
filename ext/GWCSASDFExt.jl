module GWCSASDFExt

using ASDF
using GWCS
using PrecompileTools: @compile_workload

GWCS._tag(node::ASDF.TaggedMapping) = node.tag
GWCS._tag(node::ASDF.TaggedSequence) = node.tag
GWCS._tag(node::ASDF.TaggedScalar) = node.tag
GWCS._unwrap(node::ASDF.TaggedScalar) = String(node)

function GWCS._materialize_array(array::ASDF.NDArray)
    value = array[]
    return value isa AbstractMatrix ? permutedims(value) : value
end

_asdf_metadata(value::GWCS.TaggedMapping) = ASDF.TaggedMapping(value.tag, Dict{Any,Any}(_asdf_metadata(key) => _asdf_metadata(item) for (key, item) in value))
_asdf_metadata(value::GWCS.TaggedSequence) = ASDF.TaggedSequence(value.tag, map(_asdf_metadata, value.value))
_asdf_metadata(value::GWCS.TaggedScalar) = ASDF.TaggedScalar(value.tag, String(value.value))
_asdf_metadata(value::AbstractDict) = Dict{Any,Any}(_asdf_metadata(key) => _asdf_metadata(item) for (key, item) in value)
_asdf_metadata(value::AbstractArray) = map(_asdf_metadata, value)
_asdf_metadata(value::Tuple) = map(_asdf_metadata, value)
_asdf_metadata(value::Union{GWCS.WCS,GWCS.Step,GWCS.AbstractTransform,GWCS.AbstractCoordinateFrame,GWCS.AbstractBoundingBox,GWCS.AbstractLabelMapper}) = _asdf_metadata(GWCS.to_metadata(value))
_asdf_metadata(value) = value

function GWCS.write_asdf(filename, wcs::GWCS.WCS; key="wcs")
    return ASDF.write_file(filename, Dict(String(key) => _asdf_metadata(wcs)))
end

function GWCS.write_asdf(filename, document::AbstractDict)
    return ASDF.write_file(filename, _asdf_metadata(document))
end

# Parses the main workload's pipeline as the OrderedDict-backed nodes `ASDF.load_file` builds, so the cached
# adapter and parser specializations match loaded files.
_ordered(value::ASDF.TaggedMapping) = ASDF.TaggedMapping(value.tag, ASDF.OrderedDict{Any,Any}(_ordered(k) => _ordered(v) for (k, v) in value.value))
_ordered(value::ASDF.TaggedSequence) = ASDF.TaggedSequence(value.tag, map(_ordered, value.value))
_ordered(value::AbstractDict) = ASDF.OrderedDict{Any,Any}(_ordered(k) => _ordered(v) for (k, v) in value)
_ordered(value::AbstractArray) = map(_ordered, value)
_ordered(value) = value
@compile_workload begin
    wcs = GWCS.parse_wcs(_ordered(_asdf_metadata(GWCS._load_tagged_yaml(GWCS._WORKLOAD_ASDF))))
    GWCS.world_to_pixel(wcs, GWCS.pixel_to_world(wcs, 100.0, 200.0))
end

end
