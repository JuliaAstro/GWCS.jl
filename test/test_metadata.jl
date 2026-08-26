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
end
