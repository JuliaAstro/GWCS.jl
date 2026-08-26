@testset "Numerical inverse" begin
    pixels = Frame2D(name="pixels")
    world = CoordinateFrame{2}(name="world")
    distortion = parallel(Polynomial1D([0.0, 1.0, 1e-3]), Polynomial1D([0.0, 1.0, -5e-4]))
    wcs = WCS(Step(pixels, distortion), Step(world); bounding_box=BoundingBox(((-10.0, 10.0), (-10.0, 10.0))))

    pixel = SVector(2.5, -3.5)
    transformed = pixel_to_world(wcs, pixel)
    @test wcs.backward isa NoInverse
    @test world_to_pixel(wcs, transformed; tolerance=1e-10, quiet=false) ≈ pixel atol=1e-9
    @test numerical_inverse(wcs, transformed; tolerance=1e-10, quiet=false) ≈ pixel atol=1e-9

    best = world_to_pixel(wcs, transformed; tolerance=1e-15, maxiter=1)
    @test all(isfinite, best)
    error = try
        world_to_pixel(wcs, transformed; tolerance=1e-15, maxiter=1, quiet=false)
        nothing
    catch exception
        exception
    end
    @test error isa NoConvergence
    @test error.niter == 1
    @test length(error.best_solution) == 2

    non_square = WCS(Step(CoordinateFrame{2}(name="input"), Polynomial2D([0.0 1.0; 1.0 0.0])), Step(CoordinateFrame{1}(name="output")))
    @test_throws "numerical inversion requires equal pixel and world dimensions" world_to_pixel(non_square, 1.0)

    permuted_pixels = CoordinateFrame{3}(name="permuted_pixels", axes_order=(2, 3, 1))
    permuted_world = CoordinateFrame{3}(name="permuted_world", axes_order=(3, 1, 2))
    permuted_distortion = parallel(Polynomial1D([0.0, 1.0, 1e-3]), Polynomial1D([0.0, 1.0, -5e-4]), Polynomial1D([0.0, 1.0, 2e-4]))
    permuted_wcs = WCS(Step(permuted_pixels, permuted_distortion), Step(permuted_world);
        bounding_box=BoundingBox(((-10.0, 10.0), (-10.0, 10.0), (-10.0, 10.0))))
    permuted_pixel = SVector(2.5, -3.5, 1.25)
    permuted_output = pixel_to_world(permuted_wcs, permuted_pixel)
    @test permuted_wcs.backward isa NoInverse
    @test world_to_pixel(permuted_wcs, permuted_output; tolerance=1e-10, quiet=false) ≈ permuted_pixel atol=1e-9

    arithmetic_distortion = IdentityTransform{1}() + Polynomial1D([0.0, 0.0, 1e-3])
    arithmetic_wcs = WCS(Step(CoordinateFrame{1}(name="arithmetic_input"), arithmetic_distortion), Step(CoordinateFrame{1}(name="arithmetic_output")); bounding_box=BoundingBox(((-10.0, 10.0),)))
    arithmetic_pixel = SVector(2.5)
    arithmetic_world = pixel_to_world(arithmetic_wcs, arithmetic_pixel)
    @test arithmetic_wcs.backward isa NoInverse
    @test world_to_pixel(arithmetic_wcs, arithmetic_world; tolerance=1e-12, quiet=false) ≈ arithmetic_pixel atol=1e-12
end
