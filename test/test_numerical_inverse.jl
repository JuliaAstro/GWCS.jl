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
end
