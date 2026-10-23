using ASDF
using Logging

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
    # in_image reads the solver's convergence status instead of catching a NoConvergence, which allocates.
    allocated_in_image(wcs, value) = @allocated in_image(wcs, value; tolerance=1e-15, maxiter=1)
    allocated_in_image(wcs, transformed)
    @test allocated_in_image(wcs, transformed) == 0
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

    # Newton runs only on the non-invertible prefix; the invertible tail is inverted once and never evaluated forward.
    tail_calls = Ref(0)
    counted_tail = CallableTransform{2,2}((x, y) -> (tail_calls[] += 1; (2x + 1, 3y - 1)); inverse=(u, v) -> ((u - 1) / 2, (v + 1) / 3))
    split_wcs = WCS(Step(pixels, distortion | counted_tail), Step(world))
    split_world = pixel_to_world(split_wcs, pixel)
    tail_calls[] = 0
    @test world_to_pixel(split_wcs, split_world; tolerance=1e-12) ≈ pixel atol=1e-12
    @test tail_calls[] == 0

    sip_header = Dict("NAXIS" => 2, "CTYPE1" => "RA---TAN-SIP", "CTYPE2" => "DEC--TAN-SIP", "CRPIX1" => 512.0, "CRPIX2" => 512.0, "CRVAL1" => 83.8, "CRVAL2" => -5.4,
        "CD1_1" => -2.7778e-4, "CD1_2" => 0.0, "CD2_1" => 0.0, "CD2_2" => 2.7778e-4, "A_ORDER" => 2, "A_2_0" => 5e-6, "B_ORDER" => 2, "B_0_2" => 2e-6)
    sip = parse_fits_wcs(sip_header)
    sip_world = pixel_to_world(sip, SVector(100.0, 900.0))
    @test @inferred(world_to_pixel(sip, sip_world)) ≈ SVector(100.0, 900.0) atol=1e-9
    allocated_inverse(wcs, value) = @allocated world_to_pixel(wcs, value)
    allocated_inverse(sip, sip_world)
    @test allocated_inverse(sip, sip_world) == 0

    # With an image shape, construction fits pixels as a polynomial of the prefix output, which starts Newton within
    # 0.01 px of an in-image root; without one, or with an AP/BP reverse, there is no fit.
    has_fit(wcs) = all(isfinite, wcs.inverse_split.approximate.low)
    @test !has_fit(sip)
    shaped = parse_fits_wcs(merge(sip_header, Dict("NAXIS1" => 1024, "NAXIS2" => 1024)))
    @test has_fit(shaped)
    shaped_intermediate = shaped.inverse_split.tail_inverse(sip_world)
    @test maximum(abs, GWCS._initial_guess(shaped, shaped_intermediate, Float64, nothing, shaped.inverse_split.prefix) - SVector(100.0, 900.0)) < 0.01
    @test GWCS._numerical_result(shaped, sip_world)[4] == 1
    @test @inferred(world_to_pixel(shaped, sip_world)) ≈ SVector(100.0, 900.0) atol=1e-9
    allocated_inverse(shaped, sip_world)
    @test allocated_inverse(shaped, sip_world) == 0
    @test !has_fit(parse_fits_wcs(merge(sip_header, Dict("NAXIS1" => 1024, "NAXIS2" => 1024, "AP_ORDER" => 2, "AP_2_0" => -5e-6, "BP_ORDER" => 2, "BP_0_2" => -2e-6))))
    # The fitted degree, and whether a fit exists at all, depend on the coefficients, so they must not change the WCS
    # type: otherwise every file would compile its own solver, which the precompile workload cannot cover.
    stronger = parse_fits_wcs(merge(sip_header, Dict("NAXIS1" => 1024, "NAXIS2" => 1024, "A_2_0" => 2e-4)))
    folded = parse_fits_wcs(merge(sip_header, Dict("NAXIS1" => 1024, "NAXIS2" => 1024, "A_2_0" => 2e-3)))
    @test degree(stronger.inverse_split.approximate.fit.x) > degree(shaped.inverse_split.approximate.fit.x)
    @test !has_fit(folded)
    @test typeof(stronger) == typeof(shaped) == typeof(folded)
    # Outside the image the fit would extrapolate, so the usual guess applies.
    far_world = pixel_to_world(sip, SVector(6000.0, -4000.0))
    far_intermediate = shaped.inverse_split.tail_inverse(far_world)
    @test GWCS._initial_guess(shaped, far_intermediate, Float64, nothing, shaped.inverse_split.prefix) == GWCS._initial_guess(sip, far_intermediate, Float64, nothing, sip.inverse_split.prefix)
    @test world_to_pixel(shaped, far_world; with_bounding_box=false) ≈ world_to_pixel(sip, far_world; with_bounding_box=false) atol=1e-9

    # Invertible prefix leaves, such as a selector with an inverse label mapper, carry precomputed inverses for the
    # initial guess, which would otherwise build an inverse selector on every solve.
    paired = WithInverseLabelMapper(CallableLabelMapper{2,Int}((x, y) -> x < 5 ? 1 : 2), CallableLabelMapper{2,Int}((x, y) -> x < 6 ? 1 : 2))
    regions = RegionSelector(paired, (1 => (Shift(1.0) & Shift(1.0)), 2 => (Shift(2.0) & Shift(2.0))))
    selector_wcs = WCS(Step(pixels, regions | distortion), Step(world))
    selector_world = pixel_to_world(selector_wcs, pixel)
    @test selector_wcs.backward isa NoInverse
    @test @inferred(world_to_pixel(selector_wcs, selector_world)) ≈ pixel atol=1e-9
    allocated_inverse(selector_wcs, selector_world)
    @test allocated_inverse(selector_wcs, selector_world) == 0

    # Callable transforms turn a DomainError from user code (here Base sqrt) into NaN, like GWCS transforms.
    @test all(isnan, CallableTransform{1,1}(x -> sqrt(x))(SVector(-1.0)))
    # The first Newton step from 1 lands at -0.8, outside sqrt's domain; it must backtrack rather than throw.
    root_wcs = WCS(Step(CoordinateFrame{1}(name="root_input"), CallableTransform{1,1}(x -> sqrt(x))), Step(CoordinateFrame{1}(name="root_output")))
    @test world_to_pixel(root_wcs, SVector(0.1); initial_guess=(1.0,), tolerance=1e-12) ≈ SVector(0.01) atol=1e-12
end

@testset "Numerical inverse world-space fallback" begin
    # A Python GWCS pipeline (validation/generate_sky_correction_fixture.py) whose only non-invertible step is a
    # sky-space correction after the projection. Near RA 0/360 the corrected longitude can sit 360 degrees from the
    # inverted tail's, so the intermediate-space solve fails and the world-space solve must take over.
    fixture = with_logger(NullLogger()) do
        ASDF.load(joinpath(@__DIR__, "fixtures", "sky_correction.asdf"); extensions=true)
    end
    wcs = parse_wcs(fixture["wcs"])
    pixels, world = fixture["pixels"][], fixture["world"][]
    @test wcs.backward isa NoInverse
    @test maximum(abs, pixel_to_world(wcs, pixels) - world) < 1e-12
    @test maximum(abs, world_to_pixel(wcs, world; tolerance=1e-10) - pixels) < 1e-8

    seam = findfirst(i -> pixels[:, i] == [512.0, 512.0], axes(pixels, 2))
    target = SVector(world[1, seam], world[2, seam])
    intermediate = wcs.inverse_split.tail_inverse(target)
    guess = GWCS._initial_guess(wcs, intermediate, Float64, nothing, wcs.inverse_split.prefix)
    @test !first(GWCS._newton(wcs.inverse_split.prefix, intermediate, (false, false), guess, 1e-10, 30, 12))
    @test world_to_pixel(wcs, target; tolerance=1e-10) ≈ SVector(512.0, 512.0) atol=1e-8
end
