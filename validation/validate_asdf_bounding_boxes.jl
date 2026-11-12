using ASDF
using GWCS
using GWCS: BoundingBox, CompoundBoundingBox, bounded_axes, ignored_axes, selector_axes
using StaticArrays

path = joinpath(@__DIR__, "..", "test", "fixtures", "bounding_boxes.asdf")
asdf = ASDF.load(path; extensions=true, validate_checksum=false)
simple = parse_wcs(asdf["simple"])
compound = parse_wcs(asdf["compound"])

@assert simple.bounding_box isa BoundingBox{3}
@assert bounded_axes(simple.bounding_box) == (1, 2)
@assert ignored_axes(simple.bounding_box) == (3,)
@assert pixel_to_world(simple, 2.0, 12.0, 100.0) == SVector(2.0, 12.0, 100.0)
@assert all(isnan, pixel_to_world(simple, 7.0, 12.0, 100.0))

@assert compound.bounding_box isa CompoundBoundingBox{3}
@assert selector_axes(compound.bounding_box) == (3,)
@assert pixel_to_world(compound, 2.0, 3.0, 1.0) == SVector(2.0, 3.0, 1.0)
@assert all(isnan, pixel_to_world(compound, 2.0, 3.0, 2.0))

println("ASDF simple and compound bounding-box validation passed")
