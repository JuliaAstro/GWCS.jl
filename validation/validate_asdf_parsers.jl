using ASDF
using GWCS
using StaticArrays

path = joinpath(@__DIR__, "..", "test", "fixtures", "affine.asdf")
asdf = ASDF.load(path; extensions=true, validate_checksum=false)
transform = parse_transform(asdf["transform"])

@assert transform isa AnnotatedTransform
@assert transform.transform.matrix == @SMatrix [1.0 0.2; -0.1 2.0]
@assert transform.transform.offset == SVector(3.0, 4.0)
@assert transform(SVector(1.0, 2.0)) ≈ SVector(4.4, 7.9)

println("ASDF affine parser validation passed")
