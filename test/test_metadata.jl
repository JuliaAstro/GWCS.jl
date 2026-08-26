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
