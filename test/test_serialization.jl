using ASDF
using LinearAlgebra: I
using Logging
using YAML

abstract type AbstractSerializationTransform <: AbstractTransform{1,1} end

struct SerializationTransform <: AbstractSerializationTransform
    offset::Float64
end

@inline (transform::SerializationTransform)(coords::StaticVector{1}) = SVector(coords[1] + transform.offset)

function serialization_transform_parser(node, context)
    return SerializationTransform(Float64(metadata_value(required_metadata(node, "offset", context))))
end

@testset "Serialization registry" begin
    registry = TransformRegistry()
    @test_throws "no metadata serializer is registered" to_metadata(Shift(1.0); registry)
    tag = "tag:example.org:gwcs/serialization_transform"
    register_transform!(registry, tag, serialization_transform_parser; versions=(v"1.0.0",))
    register_transform_serializer!(registry, AbstractSerializationTransform, tag, (transform, context) -> Dict("offset" => -1.0); version=v"1.0.0")
    register_transform_serializer!(registry, SerializationTransform, tag, (transform, context) -> Dict("offset" => transform.offset); version=v"1.0.0")

    transform = SerializationTransform(3.0)
    node = to_metadata(transform; registry)
    @test node.tag == "$tag-1.0.0"
    @test node["offset"] == 3.0
    @test parse_transform(node; registry)(SVector(2.0)) == SVector(5.0)

    copied = copy(registry)
    @test to_metadata(transform; registry=copied)["offset"] == 3.0
    @test_throws "already registered" register_transform_serializer!(registry, SerializationTransform, tag, (transform, context) -> Dict(); version=v"1.0.0")
    @test_throws "has no matching parser registration" register_transform_serializer!(registry, Shift, "tag:example.org:gwcs/missing", (transform, context) -> Dict(); version=v"1.0.0")
    @test_throws "no metadata serializer is registered" to_metadata(CallableTransform(x -> x, 1, 1))
end

function test_transform_roundtrip(transform, coords; atol=0, rtol=0)
    node = to_metadata(transform)
    parsed = parse_transform(node)
    @test parsed(coords) ≈ transform(coords) atol=atol rtol=rtol
    return parsed, node
end

@testset "Transform serialization" begin
    test_transform_roundtrip(IdentityTransform{2}(), SVector(1.0, 2.0))
    test_transform_roundtrip(Shift(2.0), SVector(3.0))
    test_transform_roundtrip(Scale(2.0), SVector(3.0))
    test_transform_roundtrip(AffineTransform([1.0 0.2; -0.1 2.0], [3.0, 4.0]), SVector(1.0, 2.0))
    test_transform_roundtrip(AffineTransform(reshape([2.0], 1, 1), [3.0]), SVector(4.0))
    test_transform_roundtrip(ConstantTransform(4.0), SVector{0,Float64}())
    test_transform_roundtrip(ConstantTransform(4.0; n_inputs=2), SVector(1.0, 2.0))
    test_transform_roundtrip(RemapAxes{3}((3, 1)), SVector(1.0, 2.0, 3.0))
    test_transform_roundtrip(Polynomial1D([1.0, 2.0, 3.0]; domain=(-1.0, 1.0), window=(0.0, 2.0)), SVector(0.25))
    test_transform_roundtrip(Polynomial2D([1.0 2.0; 3.0 0.0]; domain=((-1.0, 1.0), (-2.0, 2.0)), window=((0.0, 1.0), (0.0, 1.0))), SVector(0.25, 0.5))
    test_transform_roundtrip(TangentPlaneToNative(), SVector(1.0, 2.0))
    test_transform_roundtrip(NativeToTangentPlane(), SVector(10.0, 45.0); atol=1e-14)
    test_transform_roundtrip(NativeToCelestial(1.0, 2.0, 3.0), SVector(10.0, 45.0); atol=1e-12)
    test_transform_roundtrip(CelestialToNative(1.0, 2.0, 3.0), SVector(10.0, 45.0); atol=1e-12)
    test_transform_roundtrip(SphericalToCartesian(180), SVector(10.0, 45.0); atol=1e-15)
    test_transform_roundtrip(CartesianToSpherical(180), SVector(1.0, 2.0, 3.0); atol=1e-15)

    rotation = Rotation3D([10.0, 20.0, 30.0], "xyz")
    parsed_rotation, rotation_node = test_transform_roundtrip(rotation, SVector(1.0, 2.0, 3.0); atol=1e-14)
    @test rotation_node.tag == "tag:stsci.edu:asdf/transform/rotate_sequence_3d-1.2.0"
    @test parsed_rotation.matrix ≈ rotation.matrix atol=1e-14

    locked = Rotation3D([10.0, 90.0, 30.0], "xyz")
    parsed_locked, _ = test_transform_roundtrip(locked, SVector(1.0, 2.0, 3.0); atol=1e-14)
    @test parsed_locked.matrix ≈ locked.matrix atol=1e-14
    near_locked = Rotation3D([10.0, 89.999, 30.0], "xyz")
    parsed_near_locked, _ = test_transform_roundtrip(near_locked, SVector(1.0, 2.0, 3.0); atol=1e-12)
    @test parsed_near_locked.matrix ≈ near_locked.matrix atol=1e-12
    @test_throws "must be orthogonal" to_metadata(Rotation3D([1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 1.0]))

    composed = SphericalToCartesian() | Rotation3D([1.0, 2.0, 3.0], "xyz") | CartesianToSpherical()
    _, composed_node = test_transform_roundtrip(composed, SVector(10.0, 45.0); atol=1e-14)
    @test length(composed_node["forward"]) == 2
    @test composed_node["forward"][2] isa GWCS.TaggedMapping
    parallel_transform = Polynomial1D([1.0, 2.0]) & Polynomial1D([3.0, 4.0]) & Polynomial1D([5.0, 6.0])
    _, parallel_node = test_transform_roundtrip(parallel_transform, SVector(2.0, 5.0, 7.0))
    @test length(parallel_node["forward"]) == 2
    @test parallel_node["forward"][2] isa GWCS.TaggedMapping
    test_transform_roundtrip(elementwise(+, Shift(1.0), Scale(2.0)), SVector(3.0))

    annotated = annotate_transform(Shift(1.0); name="offset", inputs=("pixel",), outputs=("world",))
    parsed_annotated, _ = test_transform_roundtrip(annotated, SVector(3.0))
    @test transform_name(parsed_annotated) == "offset"
    @test input_names(parsed_annotated) == ("pixel",)
    @test output_names(parsed_annotated) == ("world",)

    explicit = WithInverse(Polynomial1D([0.0, 2.0]), Scale(0.5))
    parsed_explicit, _ = test_transform_roundtrip(explicit, SVector(3.0))
    @test inverse_transform(parsed_explicit)(SVector(6.0)) == SVector(3.0)

    fixed = fix_inputs(annotate_transform(Shift(1.0) & Scale(2.0); inputs=("x", "y")), "y" => 5.0)
    parsed_fixed, _ = test_transform_roundtrip(fixed, SVector(2.0))
    @test inverse_transform(parsed_fixed)(fixed(SVector(2.0))) ≈ SVector(2.0)
    @test_throws "fixed-input atol cannot be represented" to_metadata(fix_inputs(Shift(1.0) & Scale(2.0), 2 => 5.0; atol=1.0))

    @test_throws "no metadata serializer is registered" to_metadata(AffineTransform(Matrix{Float64}(I, 3, 3)))
    @test_throws "no metadata serializer is registered" to_metadata(ConstantTransform(1.0; n_inputs=3))
end

@testset "Selector serialization" begin
    fixture = GWCS._load_tagged_yaml(read(joinpath(@__DIR__, "fixtures", "asdf_selectors.yaml"), String))
    coordinates = Dict("array" => SVector(0.0, 1.0), "dictionary" => SVector(-2.0005, 8.0), "range" => SVector(1.0, 8.0), "transform" => SVector(99.0, 8.0))
    for name in ("array", "dictionary", "range", "transform")
        original = parse_transform(fixture[name])
        parsed = parse_transform(to_metadata(original))
        @test parsed(coordinates[name]) ≈ original(coordinates[name])
    end

    mapper = ArrayLabelMapper([9 1; 2 9]; no_label=9)
    parsed_mapper = GWCS._parse_label_mapper(to_metadata(mapper), "mapper", default_transform_registry())
    @test map_label(parsed_mapper, SVector(0.0, 0.0)) === nothing
    @test map_label(parsed_mapper, SVector(1.0, 0.0)) == 1

    dictionary = DictionaryLabelMapper{1,Float64}((1.0 => ConstantTransform(2.0; n_inputs=1),); rtol=1e-4)
    @test_throws "cannot be represented" to_metadata(dictionary)
    @test_throws "CallableLabelMapper cannot be serialized" to_metadata(CallableLabelMapper{1,Int}(x -> 1))
    @test_throws "WithInverseLabelMapper cannot be represented" to_metadata(WithInverseLabelMapper(mapper, mapper))
end

@testset "Frame, domain, and WCS serialization" begin
    mixed = parse_wcs(read(joinpath(@__DIR__, "fixtures", "mixed_frame_wcs.yaml"), String))
    mixed_node = to_metadata(mixed)
    reparsed_mixed = parse_wcs(mixed_node)
    pixels = SVector(1.0, 2.0, 3.0, 4.0, 5.0)
    @test available_frames(reparsed_mixed) == available_frames(mixed)
    @test output_frame(reparsed_mixed) isa CompositeFrame
    @test output_frame(reparsed_mixed).frames[1] isa CelestialFrame
    @test output_frame(reparsed_mixed).frames[2] isa SpectralFrame
    @test output_frame(reparsed_mixed).frames[3] isa TemporalFrame
    @test output_frame(reparsed_mixed).frames[4] isa StokesFrame
    @test pixel_to_world(reparsed_mixed, pixels) == pixel_to_world(mixed, pixels)
    celestial_metadata = to_metadata(output_frame(parse_wcs(read(joinpath(@__DIR__, "fixtures", "roman_wcs.yaml"), String))))
    @test celestial_metadata["reference_frame"]["frame_attributes"] == Dict()

    bounding_fixture = GWCS._load_tagged_yaml(read(joinpath(@__DIR__, "fixtures", "asdf_bounding_boxes.yaml"), String))
    for name in ("simple_ignored", "compound")
        original = parse_wcs(bounding_fixture[name])
        reparsed = parse_wcs(to_metadata(original))
        @test typeof(reparsed.bounding_box).name.wrapper == typeof(original.bounding_box).name.wrapper
        @test in_domain(reparsed, name == "compound" ? (2.0, 3.0, 1.0) : (12.0, 2.0, 100.0))
    end

    roman = parse_wcs(read(joinpath(@__DIR__, "fixtures", "roman_wcs.yaml"), String))
    roman_node = to_metadata(roman)
    reparsed_roman = parse_wcs(roman_node)
    world = pixel_to_world(roman, 11.0, 22.0)
    @test roman_node.tag == "tag:stsci.edu:gwcs/wcs-1.4.0"
    @test reparsed_roman.name == roman.name
    @test reparsed_roman.pixel_shape == roman.pixel_shape
    @test reparsed_roman.bounding_box.intervals == roman.bounding_box.intervals
    @test pixel_to_world(reparsed_roman, 11.0, 22.0) ≈ world atol=1e-12
    @test world_to_pixel(reparsed_roman, world) ≈ SVector(11.0, 22.0) atol=1e-10

    yaml = YAML.write(roman_node)
    @test occursin("tag:stsci.edu:gwcs/wcs-1.4.0", yaml)

    simple = WCS(Step(Frame2D(; name="detector"), AffineTransform([2.0 0.0; 0.0 3.0], [1.0, -1.0])), Step(Frame2D(; name="world")))
    reparsed_simple = parse_wcs(to_metadata(simple))
    serialized_forward(wcs::W, coords::C) where {W<:WCS,C<:StaticVector} = pixel_to_world(wcs, coords)
    allocated_serialized_forward(wcs::W, coords::C) where {W<:WCS,C<:StaticVector} = @allocated pixel_to_world(wcs, coords)
    coords = SVector(2.0, 4.0)
    pixel_to_world(reparsed_simple, coords)
    @test @inferred(serialized_forward(reparsed_simple, coords)) == SVector(5.0, 11.0)
    @test allocated_serialized_forward(reparsed_simple, coords) == 0
end

@testset "ASDF writing" begin
    original = parse_wcs(read(joinpath(@__DIR__, "fixtures", "roman_wcs.yaml"), String))
    mktempdir() do directory
        filename = joinpath(directory, "wcs.asdf")
        write_asdf(filename, original)
        loaded = with_logger(NullLogger()) do
            ASDF.load(filename; extensions=true)
        end
        reparsed = parse_wcs(loaded["wcs"])
        @test pixel_to_world(reparsed, 11.0, 22.0) ≈ pixel_to_world(original, 11.0, 22.0) atol=1e-12

        document_filename = joinpath(directory, "document.asdf")
        write_asdf(document_filename, Dict("model" => Dict("wcs" => original), "name" => "canonical"))
        document = with_logger(NullLogger()) do
            ASDF.load(document_filename; extensions=true)
        end
        @test document["name"] == "canonical"
        @test pixel_to_world(parse_wcs(document["model"]["wcs"]), 11.0, 22.0) ≈ pixel_to_world(original, 11.0, 22.0) atol=1e-12
    end
end
