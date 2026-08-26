@testset "WCS pipeline" begin
    detector = Frame2D(name="detector", axes_names=("x", "y"), unit=("pixel", "pixel"))
    focal = Frame2D(name="focal", axes_names=("u", "v"), unit=("deg", "deg"))
    world = CelestialFrame(name="world")

    detector_to_focal = parallel(Shift(-10.0), Shift(-20.0))
    focal_to_world = compose(parallel(Scale(0.1), Scale(0.1)), SphericalToCartesian(), Rotation3D([10.0], "z"), CartesianToSpherical())
    wcs = WCS(Step(detector, detector_to_focal), Step(focal, focal_to_world), Step(world); name="test", pixel_shape=(100, 200), bounding_box=BoundingBox(((-0.5, 99.5), (-0.5, 199.5))))

    @test available_frames(wcs) == ("detector", "focal", "world")
    @test input_frame(wcs) === detector
    @test output_frame(wcs) === world
    @test wcs.name == "test"
    @test wcs.pixel_shape == (100, 200)
    @test SVector(10.0, 20.0) in wcs.bounding_box
    @test !(SVector(101.0, 20.0) in wcs.bounding_box)
    c_order_box = BoundingBox(((10.0, 20.0), (1.0, 2.0)); order=:C)
    @test c_order_box.order === :F
    @test SVector(1.5, 15.0) in c_order_box
    @test !(SVector(15.0, 1.5) in c_order_box)

    focal_position = transform(wcs, "detector", "focal", 11.0, 22.0)
    @test focal_position == SVector(1.0, 2.0)
    world_position = transform(wcs, "detector", "world", 11.0, 22.0)
    @test transform(wcs, "world", "detector", world_position) ≈ SVector(11.0, 22.0) atol=1e-10
    @test get_transform(wcs, "focal", "focal") isa IdentityTransform{2}

    @test_throws "frame names must be unique" WCS(Step(detector, detector_to_focal), Step(detector))
    @test_throws DimensionMismatch WCS(Step(detector, Shift(1.0)), Step(world))

    sliced = fix_inputs(wcs, "y" => 20.0)
    @test pixel_n_dim(sliced) == 1
    @test world_n_dim(sliced) == 2
    @test input_frame(sliced).axes_names == ("x",)
    @test input_frame(sliced).unit == ("pixel",)
    @test sliced.pixel_shape == (100,)
    @test sliced.bounding_box.intervals == ((-0.5, 99.5),)
    @test pixel_to_world(sliced, 11.0) ≈ pixel_to_world(wcs, 11.0, 20.0)
    @test world_to_pixel(sliced, pixel_to_world(sliced, 11.0)) ≈ SVector(11.0) atol=1e-10
    @test_throws DomainError world_to_pixel(sliced, pixel_to_world(wcs, 11.0, 21.0))
    @test_throws "lies outside the WCS bounding box" fix_inputs(wcs, "y" => 300.0)

    constant_wcs = fix_inputs(wcs, ("x" => 11.0, "y" => 20.0))
    @test pixel_n_dim(constant_wcs) == 0
    @test input_frame(constant_wcs).axes_names == ()
    @test constant_wcs.pixel_shape == ()
    @test constant_wcs.bounding_box === nothing
    @test pixel_to_world(constant_wcs) ≈ pixel_to_world(wcs, 11.0, 20.0)
    @test pixel_to_world(constant_wcs, ()) ≈ pixel_to_world(wcs, 11.0, 20.0)
    @test world_to_pixel(constant_wcs, pixel_to_world(constant_wcs)) == SVector{0,Float64}()
    @test eltype(world_to_pixel(constant_wcs, pixel_to_world(constant_wcs))) === Float64
end
