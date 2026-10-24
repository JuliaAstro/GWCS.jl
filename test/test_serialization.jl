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
    test_transform_roundtrip(IntermediateToNative(GnomonicProjection()), SVector(1.0, 2.0))
    test_transform_roundtrip(NativeToIntermediate(GnomonicProjection()), SVector(10.0, 45.0); atol=1e-14)
    for projection in (
        ZenithalPerspective(1.5, 15.0),
        SlantZenithalPerspective(1.5, 15.0, 80.0),
        StereographicProjection(),
        SlantOrthographic(0.1, 0.2),
        ZenithalEquidistant(),
        ZenithalEqualArea(),
        AiryProjection(75.8)
    )
        test_transform_roundtrip(IntermediateToNative(projection), SVector(1.0, 2.0); atol=1e-12)
        test_transform_roundtrip(NativeToIntermediate(projection), SVector(30.0, 60.0); atol=1e-12)
    end
    azp_node = to_metadata(IntermediateToNative(ZenithalPerspective(1.5, 15.0)))
    @test azp_node.tag == "tag:stsci.edu:asdf/transform/zenithal_perspective-1.5.0"
    @test azp_node["direction"] == "pix2sky"
    @test azp_node["mu"] == 1.5
    @test azp_node["gamma"] == 15.0
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

    # parallel() fuses these into block-diagonal affines, which astropy can only represent block by block.
    _, shifts_node = test_transform_roundtrip(Shift(1.0) & Shift(2.0) & Shift(3.0), SVector(1.0, 2.0, 3.0))
    @test shifts_node.tag == "tag:stsci.edu:asdf/transform/concatenate-1.4.0"
    @test shifts_node["forward"][1].tag == "tag:stsci.edu:asdf/transform/shift-1.4.0"
    test_transform_roundtrip(AffineTransform([1.0 0.2; -0.1 2.0], [3.0, 4.0]) & Scale(2.0), SVector(1.0, 2.0, 3.0))
    test_transform_roundtrip(AffineTransform(Matrix{Float64}(I, 3, 3)), SVector(1.0, 2.0, 3.0))
    @test_throws "block diagonal with 1- and 2-dimensional blocks" to_metadata(AffineTransform([1.0 0.0 1.0; 0.0 1.0 0.0; 0.0 0.0 1.0]))
    pair = compose(RemapAxes{2}((1, 2, 1, 2)), parallel(Polynomial2D([0.1 0.2 0.3; 1.0 0.01 0.0; 0.02 0.0 0.0]), Polynomial2D([0.5 1.1 0.03; -0.2 0.0 0.0; 0.04 0.0 0.0])))
    parsed_pair, pair_node = test_transform_roundtrip(pair, SVector(0.25, 0.5))
    @test pair_node.tag == "tag:stsci.edu:asdf/transform/compose-1.4.0"
    @test GWCS._numeric_transform(parsed_pair) == pair
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
    # The stored undefined value round-trips exactly even though evaluation returns NaN.
    @test GWCS.metadata_value(to_metadata(parse_transform(fixture["array"])))["undefined_transform_value"] == -100

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
    # An astropy frame other than a bare ICRS is stored as parsed, so it is written back with its tags.
    roman = read(joinpath(@__DIR__, "fixtures", "roman_wcs.yaml"), String)
    bare_icrs = "!<tag:astropy.org:astropy/coordinates/frames/icrs-1.1.0>\n        frame_attributes: {}"
    stored(reference) = to_metadata(output_frame(parse_wcs(replace(roman, bare_icrs => reference))))["reference_frame"]
    galactic = stored("!<tag:astropy.org:astropy/coordinates/frames/galactic-1.2.0>\n        frame_attributes: {}")
    @test GWCS._tag(galactic) == "tag:astropy.org:astropy/coordinates/frames/galactic-1.2.0"
    @test galactic["frame_attributes"] == Dict()
    fiducial = stored("!<tag:astropy.org:astropy/coordinates/frames/icrs-1.1.0>\n        frame_attributes: {}\n        data: !<tag:astropy.org:astropy/coordinates/representation-1.2.0> {type: UnitSphericalRepresentation, components: {lat: !<tag:astropy.org:astropy/coordinates/latitude-1.2.0> {unit: !unit/unit-1.0.0 deg, value: 45.0}}}")
    @test GWCS._tag(fiducial) == "tag:astropy.org:astropy/coordinates/frames/icrs-1.1.0"
    @test GWCS._tag(fiducial["data"]["components"]["lat"]) == "tag:astropy.org:astropy/coordinates/latitude-1.2.0"
    @test GWCS._tag(fiducial["data"]["components"]["lat"]["unit"]) == "tag:stsci.edu:asdf/unit/unit-1.0.0"
    @test fiducial["data"]["components"]["lat"]["value"] == 45.0

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
    partial_shape = with_pixel_shape(roman, (20, nothing))
    @test !haskey(to_metadata(partial_shape), "pixel_shape")
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
    mktempdir(@__DIR__) do directory
        filename = joinpath(directory, "wcs.asdf")
        write_asdf(filename, original)
        loaded = with_logger(NullLogger()) do
            ASDF.load(filename; extensions=true)
        end
        reparsed = parse_wcs(loaded["wcs"])
        @test pixel_to_world(reparsed, 11.0, 22.0) ≈ pixel_to_world(original, 11.0, 22.0) atol=1e-12
        # Annotated pipeline files keep their stored tree; only the compiled kernel fuses the linear distortion polynomials into affines.
        @test typeof(reparsed) == typeof(original)
        @test occursin("Polynomial2D", string(typeof(original.steps)))
        @test !occursin("Polynomial2D", string(typeof(original.forward)))

        document_filename = joinpath(directory, "document.asdf")
        write_asdf(document_filename, Dict("model" => Dict("wcs" => original), "name" => "canonical"))
        document = with_logger(NullLogger()) do
            ASDF.load(document_filename; extensions=true)
        end
        @test document["name"] == "canonical"
        @test pixel_to_world(parse_wcs(document["model"]["wcs"]), 11.0, 22.0) ≈ pixel_to_world(original, 11.0, 22.0) atol=1e-12
    end
end

@testset "Step frames given as names" begin
    # validation/generate_named_frames_fixture.py: Python GWCS writes these step frames as bare names (its EmptyFrame).
    fixture = with_logger(NullLogger()) do
        ASDF.load(joinpath(@__DIR__, "fixtures", "named_frames.asdf"); extensions=true)
    end
    celestial, named = parse_wcs(fixture["celestial"]), parse_wcs(fixture["named"])
    @test input_frame(celestial) isa EmptyFrame{2} && input_frame(celestial).name == "detector"
    @test output_frame(celestial) isa CelestialFrame
    @test input_frame(named) isa EmptyFrame{2} && output_frame(named) isa EmptyFrame{2} && output_frame(named).name == "focal"
    @test axis_names(input_frame(named)) == ("", "")
    @test maximum(abs, pixel_to_world(celestial, fixture["celestial_inputs"][]) - fixture["celestial_outputs"][]) < 1e-12
    @test pixel_to_world(named, fixture["named_inputs"][]) == fixture["named_outputs"][]

    # Written back as bare names, so a GWCS.jl round trip reproduces the frames and Python reads EmptyFrames again.
    @test to_metadata(named)["steps"][1]["frame"] == "detector"
    @test to_metadata(named)["steps"][2]["frame"] == "focal"
    mktempdir(@__DIR__) do directory
        filename = joinpath(directory, "named.asdf")
        write_asdf(filename, Dict("celestial" => celestial, "named" => named))
        loaded = with_logger(NullLogger()) do
            ASDF.load(filename; extensions=true)
        end
        for (original, reparsed) in ((celestial, parse_wcs(loaded["celestial"])), (named, parse_wcs(loaded["named"])))
            @test typeof(reparsed) == typeof(original)
            @test map(step -> step.frame, reparsed.steps) == map(step -> step.frame, original.steps)
            @test pixel_to_world(reparsed, 10.0, -20.0) == pixel_to_world(original, 10.0, -20.0)
        end
    end
    @test rename_frame(named, "detector", "pixels") |> input_frame == EmptyFrame{2}("pixels")
    @test_throws "a frame given only by name needs a transform to set its dimension" parse_wcs("""
    !<tag:stsci.edu:gwcs/wcs-1.4.0>
    name: ""
    steps:
    - !<tag:stsci.edu:gwcs/step-1.3.0> {frame: detector}
    """)
end

@testset "Rotation2D" begin
    # References from astropy's Rotation2D.
    rotation = Rotation2D(30.0)
    @test rotation(SVector(1.5, -2.0)) ≈ SVector(2.299038105676658, -0.9820508075688775) atol=1e-15
    @test inverse_transform(rotation)(rotation(SVector(1.5, -2.0))) ≈ SVector(1.5, -2.0) atol=1e-15
    @test transform_jacobian(rotation, SVector(1.5, -2.0)) == rotation.matrix
    @test transform_dependencies(Rotation2D(90.0)) == SMatrix{2,2,Bool}(false, true, true, false)
    @test Rotation2D(30).angle === 30.0

    # The stored step keeps the rotation, so it is written back as rotate2d; the compiled kernel fuses it into one affine.
    chain = Shift(1.0) & Shift(2.0) | Rotation2D(-75.0) | Scale(2.0) & Scale(0.5)
    wcs = WCS(Step(Frame2D(name="detector"), chain), Step(Frame2D(name="focal")))
    @test occursin("Rotation2D", string(typeof(wcs.steps[1].transform)))
    @test wcs.forward isa AffineTransform{2,2}
    @test pixel_to_world(wcs, 3.0, 4.0) ≈ SVector(13.661662276288986, -1.1553945172705744) atol=1e-13
    @test world_to_pixel(wcs, pixel_to_world(wcs, 3.0, 4.0)) ≈ SVector(3.0, 4.0) atol=1e-13
    allocated_rotation(wcs, coords) = @allocated pixel_to_world(wcs, coords)
    allocated_rotation(wcs, SVector(3.0, 4.0))
    @test allocated_rotation(wcs, SVector(3.0, 4.0)) == 0

    mktempdir(@__DIR__) do directory
        filename = joinpath(directory, "rotation.asdf")
        write_asdf(filename, wcs)
        @test occursin("!transform/rotate2d-1.5.0", read(filename, String))
    end
    reparsed = parse_wcs(to_metadata(wcs))
    @test typeof(reparsed) == typeof(wcs)
    @test pixel_to_world(reparsed, 3.0, 4.0) == pixel_to_world(wcs, 3.0, 4.0)
    parsed = parse_transform(GWCS._load_tagged_yaml("!<tag:stsci.edu:asdf/transform/rotate2d-1.0.0> {angle: 30.0}"))
    @test parsed isa Rotation2D && parsed.angle == 30.0
    @test_throws "unsupported transform tag version tag:stsci.edu:asdf/transform/rotate2d-1.6.0" parse_transform(GWCS._load_tagged_yaml("!<tag:stsci.edu:asdf/transform/rotate2d-1.6.0> {angle: 30.0}"))
end
