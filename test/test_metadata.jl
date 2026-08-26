@testset "Metadata parsing" begin
    metadata = read(joinpath(@__DIR__, "fixtures", "roman_wcs.yaml"), String)
    wcs = parse_wcs(metadata)

    @test wcs.name == "roman_test"
    @test wcs.pixel_shape == (20, 30)
    @test available_frames(wcs) == ("detector", "v2v3", "v2v3vacorr", "world")
    @test input_frame(wcs).unit == ("pixel", "pixel")
    @test output_frame(wcs).reference_frame === :ICRS
    @test wcs.backward isa AbstractTransform

    world = pixel_to_world(wcs, 11.0, 22.0)
    @test world ≈ SVector(1.0, 2.0) atol=1e-12
    @test world_to_pixel(wcs, world) ≈ SVector(11.0, 22.0) atol=1e-10
    @test all(isnan, pixel_to_world(wcs, 21.0, 22.0))

    permuted_metadata = replace(metadata, "axes_order: [0, 1]" => "axes_order: [1, 0]"; count=1)
    permuted_wcs = parse_wcs(permuted_metadata)
    @test permuted_wcs.pixel_shape == (30, 20)
    @test permuted_wcs.bounding_box.intervals == ((-0.5, 29.5), (-0.5, 19.5))
    permuted_world = pixel_to_world(permuted_wcs, 22.0, 11.0)
    @test permuted_world ≈ world atol=1e-12
    @test world_to_pixel(permuted_wcs, permuted_world) ≈ SVector(22.0, 11.0) atol=1e-10
    @test all(isnan, pixel_to_world(permuted_wcs, 11.0, 22.0))
    permuted_slice = fix_inputs(permuted_wcs, "x" => 22.0)
    @test permuted_slice.pixel_shape == (20,)
    @test permuted_slice.bounding_box.intervals == ((-0.5, 19.5),)
    @test pixel_to_world(permuted_slice, 11.0) ≈ world atol=1e-12
    invalid_order = replace(metadata, "axes_order: [0, 1]" => "axes_order: [0, 2]"; count=1)
    @test_throws "axes_order must be a permutation of 1:2" parse_wcs(invalid_order)

    bad = replace(metadata, "transform/shift-1.4.0" => "transform/not_supported-1.4.0"; count=1)
    @test_throws "unsupported transform tag" parse_wcs(bad)

    legacy_constant = parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/constant-1.3.0>
    value: 42
    """))
    @test input_dimension(legacy_constant) == 0
    @test legacy_constant(SVector{0,Float64}()) == SVector(42)
    @test inverse_transform(legacy_constant)(SVector(42.0)) == SVector{0,Float64}()

    current_constant = parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/constant-1.6.0>
    value: 3.5
    dimensions: 2
    name: background
    inputs: [x, y]
    outputs: [value]
    """))
    @test input_dimension(current_constant) == 2
    @test current_constant(SVector(1.0, 2.0)) == SVector(3.5)
    @test transform_name(current_constant) == "background"
    @test input_names(current_constant) == ("x", "y")
    @test output_names(current_constant) == ("value",)

    fixed_metadata = """
    !<tag:stsci.edu:asdf/transform/fix_inputs-1.4.0>
    forward:
      - !<tag:stsci.edu:asdf/transform/concatenate-1.3.0>
        forward:
          - !<tag:stsci.edu:asdf/transform/shift-1.4.0> {offset: 1.0}
          - !<tag:stsci.edu:asdf/transform/shift-1.4.0> {offset: 2.0}
        name: detector_offsets
        inputs: [x, y]
        outputs: [u, v]
      - keys: [y]
        values: [10.0]
    """
    fixed = parse_transform(GWCS._load_tagged_yaml(fixed_metadata))
    @test fixed isa FixedInputsTransform{1,2,2}
    @test transform_name(fixed) == "detector_offsets"
    @test input_names(fixed) == ("x",)
    @test fixed(SVector(3.0)) == SVector(4.0, 12.0)
    @test inverse_transform(fixed)(SVector(4.0, 12.0)) == SVector(3.0)
    @test_throws DomainError inverse_transform(fixed)(SVector(4.0, 13.0))

    positional = replace(fixed_metadata, "keys: [y]" => "keys: [1]")
    @test parse_transform(GWCS._load_tagged_yaml(positional))(SVector(3.0)) == SVector(4.0, 12.0)

    concatenated_default_names = """
    !<tag:stsci.edu:asdf/transform/fix_inputs-1.4.0>
    forward:
      - !<tag:stsci.edu:asdf/transform/concatenate-1.3.0>
        forward:
          - !<tag:stsci.edu:asdf/transform/shift-1.4.0> {offset: 1.0}
          - !<tag:stsci.edu:asdf/transform/shift-1.4.0> {offset: 2.0}
      - keys: [x1]
        values: [10.0]
    """
    @test parse_transform(GWCS._load_tagged_yaml(concatenated_default_names))(SVector(3.0)) == SVector(4.0, 12.0)

    implicit_names = """
    !<tag:stsci.edu:asdf/transform/fix_inputs-1.4.0>
    forward:
      - !<tag:stsci.edu:asdf/transform/polynomial-1.2.0>
        coefficients: [[1.0, 2.0], [3.0, 0.0]]
      - keys: [x]
        values: [2.0]
    """
    implicit = parse_transform(GWCS._load_tagged_yaml(implicit_names))
    @test implicit(SVector(4.0)) == SVector(15.0)

    bad_dimensions = """
    !<tag:stsci.edu:asdf/transform/constant-1.6.0>
    value: 1
    dimensions: 3
    """
    @test_throws "constant dimensions must be 1 or 2 at transform" parse_transform(GWCS._load_tagged_yaml(bad_dimensions))
    @test_throws "keys and values must have equal lengths" parse_transform(GWCS._load_tagged_yaml(replace(fixed_metadata, "values: [10.0]" => "values: []")))
    @test_throws "is not among" parse_transform(GWCS._load_tagged_yaml(replace(fixed_metadata, "keys: [y]" => "keys: [z]")))
    duplicate = replace(fixed_metadata, "keys: [y]" => "keys: [0, x]", "values: [10.0]" => "values: [10.0, 11.0]")
    @test_throws "same input cannot be fixed more than once" parse_transform(GWCS._load_tagged_yaml(duplicate))
end

@testset "ASDF arithmetic transforms" begin
    function arithmetic_metadata(name; version="1.4.0")
        return """
        !<tag:stsci.edu:asdf/transform/$name-$version>
        forward:
          - !<tag:stsci.edu:asdf/transform/shift-1.4.0> {offset: 1.0}
          - !<tag:stsci.edu:asdf/transform/scale-1.4.0> {factor: 2.0}
        """
    end

    expected = Dict("add" => 10.0, "subtract" => -2.0, "multiply" => 24.0, "divide" => 2 / 3, "power" => 4096.0)
    for version in ("1.0.0", "1.4.0"), (name, value) in expected
        parsed = parse_transform(GWCS._load_tagged_yaml(arithmetic_metadata(name; version)))
        @test parsed isa ArithmeticTransform
        @test parsed(SVector(3.0)) == SVector(value)
    end

    annotated = parse_transform(GWCS._load_tagged_yaml(arithmetic_metadata("add") * """
    name: combined
    inputs: [detector]
    outputs: [world]
    inverse: !<tag:stsci.edu:asdf/transform/compose-1.3.0>
      forward:
        - !<tag:stsci.edu:asdf/transform/shift-1.4.0> {offset: -1.0}
        - !<tag:stsci.edu:asdf/transform/scale-1.4.0> {factor: 0.3333333333333333}
    """))
    @test transform_name(annotated) == "combined"
    @test input_names(annotated) == ("detector",)
    @test output_names(annotated) == ("world",)
    @test inverse_transform(annotated)(SVector(10.0)) ≈ SVector(3.0)

    fixture = GWCS._load_tagged_yaml(read(joinpath(@__DIR__, "fixtures", "asdf_arithmetic.yaml"), String))
    for (name, value) in expected
        @test parse_transform(fixture["arithmetic"][name])(SVector(3.0)) == SVector(value)
    end
    @test parse_transform(fixture["arithmetic"]["nested"])(SVector(3.0)) == SVector(-20.0)

    parsed = parse_transform(fixture["arithmetic"]["multiply"])
    allocated_parsed_arithmetic(transform, coordinate) = @allocated transform(coordinate)
    parsed(SVector(3.0))
    @test @inferred(parsed(SVector(3.0))) == SVector(24.0)
    @test allocated_parsed_arithmetic(parsed, SVector(3.0)) == 0

    @test_throws "expected an array containing two transforms at transform/forward" parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/add-1.4.0>
    forward: invalid
    """))
    @test_throws "must contain exactly two transforms at transform/forward" parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/add-1.4.0>
    forward:
      - !<tag:stsci.edu:asdf/transform/shift-1.4.0> {offset: 1.0}
    """))
    @test_throws "must contain exactly two transforms at transform/forward" parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/add-1.4.0>
    forward:
      - !<tag:stsci.edu:asdf/transform/shift-1.4.0> {offset: 1.0}
      - !<tag:stsci.edu:asdf/transform/shift-1.4.0> {offset: 2.0}
      - !<tag:stsci.edu:asdf/transform/shift-1.4.0> {offset: 3.0}
    """))
    @test_throws "input dimensions 1 and 2 at transform" parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/add-1.4.0>
    forward:
      - !<tag:stsci.edu:asdf/transform/shift-1.4.0> {offset: 1.0}
      - !<tag:stsci.edu:asdf/transform/identity-1.4.0> {n_dims: 2}
    """))
    for name in keys(expected)
        @test_throws "unsupported transform tag version" parse_transform(GWCS._load_tagged_yaml(arithmetic_metadata(name; version="1.5.0")))
    end
end

@testset "ASDF transform parser parity" begin
    default_identity = parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/identity-1.0.0> {}
    """))
    @test default_identity isa IdentityTransform{1}
    @test default_identity(SVector(2.0)) == SVector(2.0)

    identity_transform = parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/identity-1.4.0>
    n_dims: 2
    name: passthrough
    inputs: [x, y]
    outputs: [x, y]
    """))
    @test identity_transform isa AnnotatedTransform
    @test identity_transform(SVector(1.0, 2.0)) == SVector(1.0, 2.0)
    @test transform_name(identity_transform) == "passthrough"
    @test input_names(identity_transform) == ("x", "y")

    affine_transform = parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/affine-1.5.0>
    matrix: [[1.0, 0.2], [-0.1, 2.0]]
    translation: [3.0, 4.0]
    """))
    @test affine_transform isa AffineTransform
    @test affine_transform(SVector(1.0, 2.0)) ≈ SVector(4.4, 7.9)
    @test inverse_transform(affine_transform)(SVector(4.4, 7.9)) ≈ SVector(1.0, 2.0)

    untranslated_affine = parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/affine-1.0.0>
    matrix: [[2, 0], [0, 3]]
    """))
    @test untranslated_affine(SVector(1.0, 2.0)) == SVector(2.0, 6.0)

    pix2sky = parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/gnomonic-1.4.0> {}
    """))
    sky2pix = parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/gnomonic-1.4.0>
    direction: sky2pix
    """))
    @test pix2sky isa TangentPlaneToNative
    @test sky2pix isa NativeToTangentPlane
    @test pix2sky(SVector(0.0, 0.0)) == SVector(0.0, 90.0)
    @test sky2pix(pix2sky(SVector(1.0, 2.0))) ≈ SVector(1.0, 2.0) atol=1e-14

    native_to_celestial = parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/rotate3d-1.5.0>
    phi: 10.0
    theta: 20.0
    psi: 180.0
    direction: native2celestial
    """))
    celestial_to_native = parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/rotate3d-1.5.0>
    phi: 10.0
    theta: 20.0
    psi: 180.0
    direction: celestial2native
    """))
    @test native_to_celestial(SVector(12.0, -30.0)) ≈ SVector(145.01018272767303, -75.24513591095288) atol=1e-12
    @test celestial_to_native(SVector(359.9, 45.0)) ≈ SVector(196.21638139950392, 63.638550548244325) atol=1e-12

    euler = parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/rotate3d-1.5.0>
    phi: 1.0
    theta: 2.0
    psi: 3.0
    direction: xyz
    """))
    @test euler(SVector(12.0, -30.0)) ≈ SVector(8.212668948253533, -28.241984629262163) atol=1e-12
    @test inverse_transform(euler)(euler(SVector(12.0, -30.0))) ≈ SVector(12.0, -30.0) atol=1e-12

    cartesian_sequence = parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/rotate_sequence_3d-1.2.0>
    angles: [1.0, 2.0, 3.0]
    axes_order: zyx
    rotation_type: cartesian
    """))
    @test cartesian_sequence isa Rotation3D

    spherical_sequence = parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/rotate_sequence_3d-1.2.0>
    angles: [1.0, 2.0, 3.0]
    axes_order: zyx
    rotation_type: spherical
    """))
    @test spherical_sequence(SVector(359.9, 45.0)) ≈ SVector(2.0762391545553562, 46.9751052168514) atol=1e-12
    @test inverse_transform(spherical_sequence)(spherical_sequence(SVector(12.0, -30.0))) ≈ SVector(12.0, -30.0) atol=1e-12

    fixture = read(joinpath(@__DIR__, "fixtures", "asdf_parser_wcs.yaml"), String)
    wcs = parse_wcs(fixture)
    @test available_frames(wcs) == ("detector", "world")
    @test pixel_to_world(wcs, 10.0, 20.0) ≈ SVector(84.0630895015013, -5.2010549426502735) atol=1e-12
    @test world_to_pixel(wcs, 84.0630895015013, -5.2010549426502735) ≈ SVector(10.0, 20.0) atol=1e-10
    @test @inferred(pixel_to_world(wcs, SVector(10.0, 20.0))) ≈ SVector(84.0630895015013, -5.2010549426502735) atol=1e-12

    allocated_parsed_forward(wcs, coords) = @allocated pixel_to_world(wcs, coords)
    allocated_parsed_inverse(wcs, coords) = @allocated world_to_pixel(wcs, coords)
    pixel_to_world(wcs, SVector(10.0, 20.0))
    world_to_pixel(wcs, SVector(84.0630895015013, -5.2010549426502735))
    @test allocated_parsed_forward(wcs, SVector(10.0, 20.0)) == 0
    @test allocated_parsed_inverse(wcs, SVector(84.0630895015013, -5.2010549426502735)) == 0

    @test_throws "identity dimension must be positive at transform/n_dims" parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/identity-1.4.0> {n_dims: 0}
    """))
    @test_throws "input names have length 1, expected 2 at transform/inputs" parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/identity-1.4.0>
    n_dims: 2
    inputs: [x]
    """))
    @test_throws "affine matrix must be 2 x 2 at transform/matrix" parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/affine-1.5.0>
    matrix: [[1.0]]
    """))
    @test_throws "matrix rows must have equal lengths at transform/matrix" parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/affine-1.5.0>
    matrix: [[1.0, 0.0], [0.0]]
    """))
    @test_throws "affine translation must have length 2 at transform/translation" parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/affine-1.5.0>
    matrix: [[1.0, 0.0], [0.0, 1.0]]
    translation: [1.0]
    """))
    @test_throws "expected one of pix2sky, sky2pix, got sideways at transform/direction" parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/gnomonic-1.4.0> {direction: sideways}
    """))
    @test_throws "missing required property psi at transform" parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/rotate3d-1.5.0>
    phi: 1.0
    theta: 2.0
    direction: xyz
    """))
    @test_throws "expected a unitless numeric value at transform/phi" parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/rotate3d-1.5.0>
    phi: !<tag:stsci.edu:asdf/unit/quantity-1.2.0> {value: 1.0, unit: deg}
    theta: 2.0
    psi: 3.0
    direction: xyz
    """))
    @test_throws "rotation axes must contain only x, y, or z at transform/axes_order" parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/rotate_sequence_3d-1.2.0>
    angles: [1.0, 2.0]
    axes_order: za
    rotation_type: cartesian
    """))
    @test_throws "the number of angles must match axes_order at transform" parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/rotate_sequence_3d-1.2.0>
    angles: [1.0, 2.0]
    axes_order: xyz
    rotation_type: spherical
    """))
    @test_throws "registered ranges are [1.0.0, 1.5.0)" parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/identity-1.5.0> {n_dims: 1}
    """))
    @test_throws "registered ranges are [1.0.0, 1.6.0)" parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/rotate3d-1.6.0>
    phi: 1.0
    theta: 2.0
    psi: 3.0
    direction: xyz
    """))
    @test_throws "registered ranges are [1.0.0, 1.3.0)" parse_transform(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:asdf/transform/rotate_sequence_3d-1.3.0>
    angles: [1.0]
    axes_order: x
    rotation_type: cartesian
    """))
end

@testset "Frame metadata parsing" begin
    metadata = read(joinpath(@__DIR__, "fixtures", "mixed_frame_wcs.yaml"), String)
    wcs = parse_wcs(metadata)
    output = output_frame(wcs)

    @test input_frame(wcs) isa CoordinateFrame{5}
    @test pixel_axis_types(wcs) == ntuple(_ -> PIXEL, 5)
    @test output isa CompositeFrame{5}
    @test output.axes_order == (2, 5, 1, 4, 3)
    @test world_axis_names(wcs) == ("lon", "lat", "wave", "time", "stokes")
    @test world_axis_types(wcs) == (SPATIAL, SPATIAL, SPECTRAL, TIME, STOKES)
    @test world_axis_physical_types(wcs) == ("pos.eq.ra", "pos.eq.dec", "em.wl", "time", "phys.polarization.stokes")
    @test world_axis_units(wcs) == ("deg", "deg", "nm", "s", "")
    @test output.frames[1].reference_frame === :ICRS
    @test output.frames[3].reference_epoch == 100.0
    @test output.frames[3].time_scale == "utc"
    @test output.frames[3].time_format == "jd"
    @test output.frames[3].location === nothing
    @test output.frames[3].reference_frame isa AbstractDict
    @test offset_to_time(output.frames[3], 86400.0) == 101.0
    @test axis_names(GWCS._subset_frame(output, (1, 3))) == ("lon", "wave")

    pixels = SVector(1.0, 2.0, 3.0, 4.0, 5.0)
    world = SVector(2.0, 5.0, 1.0, 4.0, 3.0)
    @test pixel_to_world(wcs, pixels) == world
    @test world_to_pixel(wcs, world) == pixels
    @test axis_correlation_matrix(wcs) == @SMatrix [false true false false false; false false false false true; true false false false false; false false false true false; false false true false false]

    bad_child = replace(metadata, "axes_names: [wave]" => "axes_names: [wave, extra]")
    @test_throws "axes_names has length 2, expected 1 at wcs/steps/2/frame/frames/2" parse_wcs(bad_child)
    duplicate_order = replace(metadata, "axes_order: [0]" => "axes_order: [1]"; count=1)
    @test_throws "axes_order entries must be unique at wcs/steps/2/frame" parse_wcs(duplicate_order)

    scalar_time = GWCS._parse_frame(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:gwcs/temporal_frame-1.2.0>
    name: time
    reference_frame: !<tag:stsci.edu:asdf/time/time-1.4.0> 2010-01-01T00:00:00
    """), "frame")
    @test scalar_time.reference_frame isa GWCS.TaggedScalar
    @test scalar_time.reference_epoch == "2010-01-01T00:00:00"
    @test_throws "does not support offset arithmetic" offset_to_time(scalar_time, 1.0)

    unnamed = GWCS._parse_frame(GWCS._load_tagged_yaml("""
    !<tag:stsci.edu:gwcs/frame-1.2.0>
    name: unnamed
    naxes: 1
    axes_type: [PIXEL]
    axes_names: [null]
    """), "frame")
    @test axis_names(unnamed) == ("",)
    @test axis_physical_types(unnamed) == ("custom:PIXEL",)
end
