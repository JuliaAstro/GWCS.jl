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
