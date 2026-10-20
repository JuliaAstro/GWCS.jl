struct ZeroTransform <: AbstractTransform{0,0} end

(::ZeroTransform)(::StaticVector{0}) = SVector{0,Float64}()

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
    @test in_image(wcs, transformed; tolerance=1e-10)

    best = world_to_pixel(wcs, transformed; tolerance=1e-15, maxiter=1, quiet=true)
    @test all(isfinite, best)
    @test !in_image(wcs, transformed; tolerance=1e-15, maxiter=1)
    error = try
        world_to_pixel(wcs, transformed; tolerance=1e-15, maxiter=1)
        nothing
    catch exception
        exception
    end
    @test error isa NoConvergence
    @test error.niter == 1
    @test length(error.best_solution) == 2
    @test_throws "maxiter must be positive" in_image(wcs, transformed; maxiter=0)

    outside_pixel = SVector(12.0, 0.0)
    outside_world = wcs.forward(outside_pixel)
    @test !in_image(wcs, outside_world; initial_guess=outside_pixel, tolerance=1e-10)
    @test in_image(wcs, outside_world; initial_guess=outside_pixel, tolerance=1e-10, with_bounding_box=false)

    non_square = WCS(Step(CoordinateFrame{2}(name="input"), Polynomial2D([0.0 1.0; 1.0 0.0])), Step(CoordinateFrame{1}(name="output")))
    @test_throws "numerical inversion requires equal pixel and world dimensions" world_to_pixel(non_square, 1.0)
    @test_throws "numerical inversion requires equal pixel and world dimensions" in_image(non_square, 1.0)

    identity_wcs = WCS(Step(CoordinateFrame{1}(name="identity_input"), CallableTransform{1,1}(identity)), Step(CoordinateFrame{1}(name="identity_output")))
    @test world_to_pixel(identity_wcs, 2.0; initial_guess=(5.0,), maxiter=1, tolerance=1e-12, quiet=false) == SVector(2.0)
    @test in_image(identity_wcs, 2.0; initial_guess=(5.0,), maxiter=1, tolerance=1e-12)

    zero_wcs = WCS(Step(CoordinateFrame{0}(name="zero_input"), ZeroTransform()), Step(CoordinateFrame{0}(name="zero_output")))
    @test world_to_pixel(zero_wcs) == SVector{0,Float64}()
    @test in_image(zero_wcs)

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

    rectangular_children = RemapAxes{2}((1, 2, 1)) | RemapAxes{3}((1, 2)) | distortion
    rectangular_wcs = WCS(Step(pixels, rectangular_children), Step(world); bounding_box=BoundingBox(((-10.0, 10.0), (-10.0, 10.0))))
    rectangular_world = pixel_to_world(rectangular_wcs, pixel)
    @test world_to_pixel(rectangular_wcs, rectangular_world; tolerance=1e-10, quiet=false) ≈ pixel atol=1e-9

    domain_initializer = GWCS._ApproximateInverseTransform(arithmetic_distortion, CallableTransform{1,1}(x -> throw(DomainError(x))))
    bounded_initializer = WCS(Step(CoordinateFrame{1}(name="bounded_input"), domain_initializer), Step(CoordinateFrame{1}(name="bounded_output")); bounding_box=BoundingBox(((2.0, 4.0),)))
    @test world_to_pixel(bounded_initializer, pixel_to_world(bounded_initializer, 2.5); tolerance=1e-12, quiet=false) ≈ SVector(2.5) atol=1e-12

    # Callable transforms turn a DomainError from user code (here Base sqrt) into NaN, like GWCS transforms.
    @test all(isnan, CallableTransform{1,1}(x -> sqrt(x))(SVector(-1.0)))
    # The first Newton step from 1 lands at -0.8, outside sqrt's domain; it must backtrack rather than throw.
    root_wcs = WCS(Step(CoordinateFrame{1}(name="root_input"), CallableTransform{1,1}(x -> sqrt(x))), Step(CoordinateFrame{1}(name="root_output")))
    @test world_to_pixel(root_wcs, SVector(0.1); initial_guess=(1.0,), tolerance=1e-12) ≈ SVector(0.01) atol=1e-12
end
