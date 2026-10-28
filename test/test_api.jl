struct OffsetColumnMatrix{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    parent::M
end

Base.size(matrix::OffsetColumnMatrix) = size(matrix.parent)
Base.axes(matrix::OffsetColumnMatrix) = (axes(matrix.parent, 1), 0:(size(matrix.parent, 2) - 1))
Base.getindex(matrix::OffsetColumnMatrix, row::Int, column::Int) = matrix.parent[row, column + 1]

@testset "Evaluation API" begin
    pixels = Frame2D(name="pixels")
    world = CoordinateFrame{2}(name="world")
    affine = AffineTransform(SMatrix{2,2}(2.0, 0.0, 0.0, 3.0), SVector(1.0, -1.0))
    wcs = WCS(Step(pixels, affine), Step(world); bounding_box=BoundingBox(((-0.5, 9.5), (-0.5, 19.5))))

    @test pixel_n_dim(wcs) == 2
    @test world_n_dim(wcs) == 2
    @test pixel_to_world(wcs, 2.0, 4.0) == SVector(5.0, 11.0)
    @test pixel_to_world(wcs, (2.0, 4.0)) == SVector(5.0, 11.0)
    @test pixel_to_world(wcs, [2.0, 4.0]) == SVector(5.0, 11.0)
    @test world_to_pixel(wcs, 5.0, 11.0) ≈ SVector(2.0, 4.0)
    @test eltype(pixel_to_world(wcs, SVector(2.0f0, 4.0f0))) === Float32
    @test eltype(world_to_pixel(wcs, SVector(5.0f0, 11.0f0))) === Float32
    @test all(isnan, pixel_to_world(wcs, 20.0, 4.0))
    @test pixel_to_world(wcs, 20.0, 4.0; with_bounding_box=false) == SVector(41.0, 11.0)
    @test all(isnan, world_to_pixel(wcs, 41.0, 11.0))
    @test in_domain(wcs, 2.0, 4.0)
    @test in_domain(wcs, (2.0, 4.0))
    @test in_domain(wcs, [2.0, 4.0])
    @test !in_domain(wcs, 20.0, 4.0)
    @test !in_domain(wcs, NaN, 4.0)
    @test in_image(wcs, 5.0, 11.0)
    @test in_image(wcs, (5.0, 11.0))
    @test in_image(wcs, [5.0, 11.0])
    @test !in_image(wcs, 41.0, 11.0)
    @test in_image(wcs, 41.0, 11.0; with_bounding_box=false)

    batch = [2.0 3.0; 4.0 5.0]
    world_batch = pixel_to_world(wcs, batch)
    @test world_batch == [5.0 7.0; 11.0 14.0]
    @test world_to_pixel(wcs, world_batch) ≈ batch
    @test in_domain(wcs, [2.0 20.0; 4.0 4.0]) == Bool[true, false]
    @test in_image(wcs, [5.0 41.0; 11.0 11.0]) == Bool[true, false]
    @test in_domain(wcs, OffsetColumnMatrix([2.0 20.0; 4.0 4.0])) == Bool[true, false]
    @test in_image(wcs, OffsetColumnMatrix([5.0 41.0; 11.0 11.0])) == Bool[true, false]
    @test in_domain(wcs, Matrix{Float64}(undef, 2, 0)) == BitVector()
    @test in_image(wcs, Matrix{Float64}(undef, 2, 0)) == BitVector()
    @test_throws DimensionMismatch pixel_to_world(wcs, [1.0])
    @test_throws DimensionMismatch world_to_pixel(wcs, [1.0])
    @test_throws DimensionMismatch in_domain(wcs, [1.0])
    @test_throws DimensionMismatch in_image(wcs, [1.0])

    coords = SVector(2.0, 4.0)
    allocated_pixel_to_world(wcs, coords) = @allocated pixel_to_world(wcs, coords)
    allocated_world_to_pixel(wcs, coords) = @allocated world_to_pixel(wcs, coords)
    allocated_in_domain(wcs, coords) = @allocated in_domain(wcs, coords)
    allocated_in_image(wcs, coords) = @allocated in_image(wcs, coords)
    pixel_to_world(wcs, coords)
    world_to_pixel(wcs, SVector(5.0, 11.0))
    @test @inferred(pixel_to_world(wcs, coords)) == SVector(5.0, 11.0)
    @test @inferred(world_to_pixel(wcs, SVector(5.0, 11.0))) ≈ coords
    @test @inferred(in_domain(wcs, coords))
    @test @inferred(in_image(wcs, SVector(5.0, 11.0)))
    @test allocated_pixel_to_world(wcs, coords) == 0
    @test allocated_world_to_pixel(wcs, SVector(5.0, 11.0)) == 0
    @test allocated_in_domain(wcs, coords) == 0
    @test allocated_in_image(wcs, SVector(5.0, 11.0)) == 0

    constant_wcs = fix_inputs(wcs, (1 => 2.0, 2 => 4.0))
    @test pixel_to_world(constant_wcs) == SVector(5.0, 11.0)
    @test world_to_pixel(constant_wcs, 5.0, 11.0) == SVector{0,Float64}()
    @test in_domain(constant_wcs)
    @test in_image(constant_wcs, 5.0, 11.0)
    @test !in_image(constant_wcs, 7.0, 11.0)
end

@testset "One-based pixels with origin=1" begin
    affine = AffineTransform(SMatrix{2,2}(2.0, 0.0, 0.0, 3.0), SVector(1.0, -1.0))
    wcs = WCS(Step(Frame2D(name="pixels"), affine), Step(CoordinateFrame{2}(name="world")); bounding_box=BoundingBox(((-0.5, 9.5), (-0.5, 19.5))))
    # Every pixel axis shifts by one; world coordinates and fill_value never do.
    @test pixel_to_world(wcs, 3.0, 5.0; origin=1) == pixel_to_world(wcs, 2.0, 4.0)
    @test pixel_to_world(wcs, [3.0 4.0; 5.0 6.0]; origin=1) == pixel_to_world(wcs, [2.0 3.0; 4.0 5.0])
    @test world_to_pixel(wcs, 5.0, 11.0; origin=1) ≈ SVector(3.0, 5.0)
    @test world_to_pixel(wcs, [5.0 7.0; 11.0 14.0]; origin=1) ≈ [3.0 4.0; 5.0 6.0]
    @test eltype(world_to_pixel(wcs, SVector(5.0f0, 11.0f0); origin=1)) === Float32
    # Pixel 10 is outside the zero-based box (-0.5, 9.5) but is its last pixel one-based.
    @test pixel_to_world(wcs, 10.0, 1.0; origin=1) == pixel_to_world(wcs, 9.0, 0.0)
    @test all(isnan, pixel_to_world(wcs, 0.0, 1.0; origin=1))
    @test pixel_to_world(wcs, 0.0, 1.0; origin=1, fill_value=-1.0) == SVector(-1.0, -1.0)
    @test world_to_pixel(wcs, 41.0, 11.0; origin=1, fill_value=-1.0) == SVector(-1.0, -1.0)
    @test in_domain(wcs, 10.0, 20.0; origin=1) && !in_domain(wcs, 10.0, 20.0)
    @test !in_domain(wcs, 0.0, 1.0; origin=1)
    @test in_domain(wcs, [10.0 0.0; 20.0 1.0]; origin=1) == Bool[true, false]
    @test grid_from_bounding_box(wcs; origin=1) == grid_from_bounding_box(wcs) .+ 1
    @test grid_from_bounding_box(wcs; fixed=(2 => 5.0), origin=1) == grid_from_bounding_box(wcs; fixed=(2 => 4.0)) .+ 1
    # Errors report the caller's pixel.
    @test_throws "pixel (11.0, 1.0) is outside the WCS bounding box" pixel_to_world(wcs, 11.0, 1.0; origin=1, strict=true)
    @test_throws "maps to pixel (21.0, 5.0), outside the WCS bounding box" world_to_pixel(wcs, 41.0, 11.0; origin=1, strict=true)
    @test_throws "origin must be 0 or 1, got 2" pixel_to_world(wcs, 1.0, 1.0; origin=2)
    @test_throws "origin must be 0 or 1, got 2" world_to_pixel(wcs, 1.0, 1.0; origin=2)
    @test_throws "origin must be 0 or 1, got 2" in_domain(wcs, 1.0, 1.0; origin=2)
    @test_throws "origin must be 0 or 1, got 2" grid_from_bounding_box(wcs; origin=2)

    # The numerical inverse also takes initial_guess one-based and reports a one-based NoConvergence solution.
    distortion = parallel(Polynomial1D([0.0, 1.0, 1e-3]), Polynomial1D([0.0, 1.0, -5e-4]))
    numerical = WCS(Step(Frame2D(name="pixels"), distortion), Step(CoordinateFrame{2}(name="world")); bounding_box=BoundingBox(((-10.0, 10.0), (-10.0, 10.0))))
    target = pixel_to_world(numerical, 2.5, -3.5)
    @test world_to_pixel(numerical, target; origin=1, tolerance=1e-12) ≈ SVector(3.5, -2.5) atol=1e-10
    @test numerical_inverse(numerical, target; origin=1, tolerance=1e-12) ≈ SVector(3.5, -2.5) atol=1e-10
    # One Newton step from the same physical pixel lands at the same physical pixel.
    one_step(origin, guess) = world_to_pixel(numerical, target; origin, initial_guess=guess, maxiter=1, quiet=true)
    @test one_step(1, (1.0, 1.0)) ≈ one_step(0, (0.0, 0.0)) .+ 1 rtol=1e-14
    @test one_step(1, (1.0, 1.0)) != one_step(0, (1.0, 1.0)) .+ 1
    @test in_image(numerical, target; origin=1, initial_guess=(3.5, -2.5))
    unconverged(origin, guess) = try
        world_to_pixel(numerical, target; origin, initial_guess=guess, tolerance=1e-30, maxiter=2)
    catch error
        error
    end
    @test unconverged(1, (1.0, 1.0)) isa NoConvergence
    @test unconverged(1, (1.0, 1.0)).best_solution ≈ unconverged(0, (0.0, 0.0)).best_solution .+ 1 rtol=1e-14

    # Strict undefined-pixel messages name the caller's pixel too.
    ait = parse_fits_wcs(Dict("NAXIS" => 2, "CTYPE1" => "RA---AIT", "CTYPE2" => "DEC--AIT", "CRPIX1" => 181.0, "CRPIX2" => 91.0, "CRVAL1" => 0.0, "CRVAL2" => 0.0, "CDELT1" => -1.0, "CDELT2" => 1.0))
    @test pixel_to_world(ait, 181.0, 91.0; origin=1) ≈ SVector(0.0, 0.0) atol=1e-12
    @test_throws r"pixel \(1\.0, 1\.0\) has no world coordinate" pixel_to_world(ait, 1.0, 1.0; origin=1, strict=true)

    coords = SVector(3.0, 5.0)
    allocated_pixel_to_world(wcs, coords) = @allocated pixel_to_world(wcs, coords; origin=1)
    allocated_world_to_pixel(wcs, coords) = @allocated world_to_pixel(wcs, coords; origin=1)
    allocated_in_domain(wcs, coords) = @allocated in_domain(wcs, coords; origin=1)
    allocated_pixel_to_world(wcs, coords); allocated_world_to_pixel(wcs, SVector(5.0, 11.0)); allocated_in_domain(wcs, coords)
    @test @inferred(pixel_to_world(wcs, coords; origin=1)) == SVector(5.0, 11.0)
    @test @inferred(world_to_pixel(wcs, SVector(5.0, 11.0); origin=1)) ≈ coords
    @test allocated_pixel_to_world(wcs, coords) == 0
    @test allocated_world_to_pixel(wcs, SVector(5.0, 11.0)) == 0
    @test allocated_in_domain(wcs, coords) == 0
    allocated_numerical(wcs, coords) = @allocated world_to_pixel(wcs, coords; origin=1)
    allocated_numerical(numerical, target)
    @test allocated_numerical(numerical, target) == 0
end

@testset "Coordinates without a result and strict evaluation" begin
    sky_header(code) = Dict("NAXIS" => 2, "CTYPE1" => "RA---$code", "CTYPE2" => "DEC--$code", "CRPIX1" => 181.0, "CRPIX2" => 91.0,
        "CRVAL1" => 0.0, "CRVAL2" => 0.0, "CDELT1" => -1.0, "CDELT2" => 1.0)
    ait = parse_fits_wcs(sky_header("AIT"))
    # Pixels are zero-based, so (180, 90) is the reference pixel CRPIX = (181, 91).
    grid = [0.0 180.0 359.0; 0.0 90.0 179.0]
    # By default, off-sky pixels of an all-sky projection return fill_value.
    filled = pixel_to_world(ait, grid)
    @test all(isnan, filled[:, 1]) && all(isnan, filled[:, 3])
    @test filled[:, 2] ≈ [0.0, 0.0] atol=1e-12
    @test pixel_to_world(ait, grid; fill_value=-999.0)[:, 1] == [-999.0, -999.0]
    @test all(isnan, pixel_to_world(ait, 0.0, 0.0))
    @test_throws DomainError pixel_to_world(ait, grid; strict=true)
    @test @inferred(pixel_to_world(ait, SVector(180.0, 90.0); strict=true)) == pixel_to_world(ait, SVector(180.0, 90.0))

    # The antipode of a TAN reference point has no pixel under the analytical inverse.
    tan = parse_fits_wcs(sky_header("TAN"))
    @test all(isnan, world_to_pixel(tan, 180.0, 0.0))
    @test_throws DomainError world_to_pixel(tan, 180.0, 0.0; strict=true)
    @test world_to_pixel(tan, [180.0 1.0; 0.0 2.0])[:, 2] ≈ world_to_pixel(tan, 1.0, 2.0)

    # Without AP/BP, TAN-SIP uses the numerical inverse; its analytical TAN tail still rejects the antipode.
    sip = parse_fits_wcs(merge(sky_header("TAN"), Dict("CTYPE1" => "RA---TAN-SIP", "CTYPE2" => "DEC--TAN-SIP", "A_ORDER" => 2, "B_ORDER" => 2, "A_2_0" => 1e-6)))
    @test sip.backward isa NoInverse
    @test all(isnan, world_to_pixel(sip, 180.0, 0.0))
    @test all(isnan, world_to_pixel(sip, 180.0, 0.0; quiet=true))
    @test_throws DomainError world_to_pixel(sip, 180.0, 0.0; strict=true)
    @test !in_image(sip, 180.0, 0.0)

    # Strict errors name the coordinate, the step, and the transform that is undefined there.
    @test_throws r"pixel \(0\.0, 0\.0\) has no world coordinate: IntermediateToNative\(HammerAitoffProjection\) in the step from frame \"intermediate\" to \"world\" is undefined at \(180\.0, -90\.0\)" pixel_to_world(ait, 0.0, 0.0; strict=true)
    @test_throws r"world coordinate \(180\.0, 0\.0\) has no pixel coordinate: NativeToIntermediate\(GnomonicProjection\) in the inverse of the step" world_to_pixel(tan, 180.0, 0.0; strict=true)
    @test_throws r"world coordinate \(180\.0, 0\.0\) has no pixel coordinate: NativeToIntermediate\(GnomonicProjection\)" world_to_pixel(sip, 180.0, 0.0; strict=true)

    # Strict mode also rejects bounding-box misses in both directions.
    boxed = WCS(Step(Frame2D(name="pixels"), AffineTransform(SMatrix{2,2}(2.0, 0.0, 0.0, 3.0), SVector(1.0, -1.0))), Step(CoordinateFrame{2}(name="world")); bounding_box=BoundingBox(((-0.5, 9.5), (-0.5, 19.5))))
    @test all(isnan, pixel_to_world(boxed, 20.0, 4.0))
    @test_throws "pixel (20.0, 4.0) is outside the WCS bounding box" pixel_to_world(boxed, 20.0, 4.0; strict=true)
    @test_throws "maps to pixel (20.0, 4.0), outside the WCS bounding box" world_to_pixel(boxed, 41.0, 11.0; strict=true)
    @test pixel_to_world(boxed, 20.0, 4.0; strict=true, with_bounding_box=false) == SVector(41.0, 11.0)

    # NaN input propagates as NaN in both modes; other errors and NoConvergence always propagate.
    @test all(isnan, pixel_to_world(ait, NaN, 1.0; strict=true))
    @test all(isnan, pixel_to_world(boxed, NaN, 1.0; strict=true))
    @test all(isnan, world_to_pixel(sip, NaN, 0.0; strict=true))
    failing = WCS(Step(Frame2D(name="pixels"), CallableTransform{2,2}((x, y) -> throw(ArgumentError("calibration unavailable")))), Step(CoordinateFrame{2}(name="world")))
    @test_throws "calibration unavailable" pixel_to_world(failing, 1.0, 2.0)
    distortion = WCS(Step(Frame2D(name="pixels"), parallel(Polynomial1D([0.0, 1.0, 1e-3]), Polynomial1D([0.0, 1.0, -5e-4]))), Step(CoordinateFrame{2}(name="world")))
    @test_throws NoConvergence world_to_pixel(distortion, 2.5, -3.5; tolerance=1e-15, maxiter=1)

    # Filling costs no more than a valid evaluation.
    allocated_default(wcs, coords) = @allocated pixel_to_world(wcs, coords)
    for coords in (SVector(1.0, 1.0), SVector(181.0, 91.0))
        pixel_to_world(ait, coords)
        @test allocated_default(ait, coords) == 0
        @test @inferred(pixel_to_world(ait, coords)) isa SVector{2,Float64}
    end
end

@testset "Pixel area" begin
    cd = [-2e-4 3e-5; 1e-5 2.5e-4]
    header(lon, lat, code="TAN") = Dict("NAXIS" => 2, "CTYPE1" => "$lon-$code", "CTYPE2" => "$lat-$code", "CRPIX1" => 51.0, "CRPIX2" => 41.0, "CRVAL1" => 30.0, "CRVAL2" => 40.0,
        "CD1_1" => cd[1, 1], "CD1_2" => cd[1, 2], "CD2_1" => cd[2, 1], "CD2_2" => cd[2, 2])
    # |det J| cos(lat) from central differences of pixel_to_world, with the latitude on world axis `latitude`.
    function difference_area(wcs, x, y; latitude=2, step=1e-3)
        dx = (pixel_to_world(wcs, x + step, y) - pixel_to_world(wcs, x - step, y)) / 2step
        dy = (pixel_to_world(wcs, x, y + step) - pixel_to_world(wcs, x, y - step)) / 2step
        return abs(dx[1] * dy[2] - dx[2] * dy[1]) * cosd(pixel_to_world(wcs, x, y)[latitude]) * deg2rad(1.0)^2
    end
    tan = parse_fits_wcs(header("RA--", "DEC-"))
    # Half a pixel from the reference pixel, the TAN area is |det CD| to well below the tolerance.
    @test pixel_area(tan, 50.5, 40.5) ≈ abs(cd[1, 1] * cd[2, 2] - cd[1, 2] * cd[2, 1]) * deg2rad(1.0)^2 rtol=1e-9
    @test pixel_area(tan, 3000.0, -2000.0) ≈ difference_area(tan, 3000.0, -2000.0) rtol=1e-7
    @test pixel_area(tan, (3000.0, -2000.0)) == pixel_area(tan, [3000.0, -2000.0]) == pixel_area(tan, SVector(3000.0, -2000.0))
    @test pixel_area(tan, 3001.0, -1999.0; origin=1) == pixel_area(tan, 3000.0, -2000.0)
    @test pixel_area(tan, [3000.0 50.5; -2000.0 40.5]) == [pixel_area(tan, 3000.0, -2000.0), pixel_area(tan, 50.5, 40.5)]
    @test pixel_area(tan, Matrix{Float64}(undef, 2, 0)) == Float64[]
    @test pixel_area(tan, NaN, 1.0) === NaN
    # FITS celestial axes become a CelestialFrame, in any order and coordinate system.
    @test pixel_area(parse_fits_wcs(header("GLON", "GLAT")), 3000.0, -2000.0) ≈ pixel_area(tan, 3000.0, -2000.0) rtol=1e-12
    swapped = parse_fits_wcs(header("DEC-", "RA--"))
    @test pixel_area(swapped, 3000.0, -2000.0) ≈ difference_area(swapped, 3000.0, -2000.0; latitude=1) rtol=1e-7
    # A CelestialFrame's latitude is its second axis; its units may be radians.
    roman = parse_wcs(read(joinpath(@__DIR__, "fixtures", "roman_wcs.yaml"), String))
    @test output_frame(roman) isa CelestialFrame
    @test pixel_area(roman, 12.0, 21.0) ≈ difference_area(roman, 12.0, 21.0) rtol=1e-7
    to_radians = parallel(Scale(deg2rad(1.0)), Scale(deg2rad(1.0)))
    radians = WCS(Step(input_frame(tan), get_transform(tan, "detector", "world") | to_radians), Step(CelestialFrame(name="world", unit=("rad", "rad"))))
    @test pixel_area(radians, 3000.0, -2000.0) ≈ pixel_area(tan, 3000.0, -2000.0) rtol=1e-12

    # Pixels outside the bounding box or off the sky are filled, or throw with strict=true.
    boxed = with_bounding_box(tan, BoundingBox(((-0.5, 99.5), (-0.5, 79.5))))
    @test isnan(pixel_area(boxed, 120.0, 10.0))
    @test pixel_area(boxed, 120.0, 10.0; fill_value=0.0) == 0.0
    @test pixel_area(boxed, 120.0, 10.0; with_bounding_box=false) == pixel_area(tan, 120.0, 10.0)
    @test_throws "pixel (120.0, 10.0) is outside the WCS bounding box" pixel_area(boxed, 120.0, 10.0; strict=true)
    ait = parse_fits_wcs(Dict("NAXIS" => 2, "CTYPE1" => "RA---AIT", "CTYPE2" => "DEC--AIT", "CRPIX1" => 181.0, "CRPIX2" => 91.0, "CDELT1" => -1.0, "CDELT2" => 1.0))
    @test isnan(pixel_area(ait, 0.0, 0.0))
    @test_throws r"pixel \(0\.0, 0\.0\) has no world coordinate: IntermediateToNative\(HammerAitoffProjection\)" pixel_area(ait, 0.0, 0.0; strict=true)
    @test pixel_area(ait, 180.5, 90.5) ≈ deg2rad(1.0)^2 rtol=1e-4

    affine = WCS(Step(Frame2D(name="pixels"), AffineTransform(SMatrix{2,2}(2.0, 0.0, 0.0, 3.0), SVector(1.0, -1.0))), Step(CoordinateFrame{2}(name="world")))
    @test_throws "pixel_area requires 2 pixel axes and a 2D CelestialFrame output frame, got 2 pixel axes and a CoordinateFrame{2}" pixel_area(affine, 1.0, 2.0)
    arcsec = WCS(Step(input_frame(tan), get_transform(tan, "detector", "world")), Step(CelestialFrame(name="world", unit=("arcsec", "arcsec"))))
    @test_throws "pixel_area requires celestial axes in deg or rad, got \"arcsec\"" pixel_area(arcsec, 1.0, 2.0)
    @test_throws DimensionMismatch pixel_area(tan, [1.0])

    allocated_pixel_area(wcs, coords) = @allocated pixel_area(wcs, coords)
    for (wcs, coords) in ((tan, SVector(3000.0, -2000.0)), (roman, SVector(12.0, 21.0)))
        @test @inferred(pixel_area(wcs, coords)) isa Float64
        @test allocated_pixel_area(wcs, coords) == 0
    end
end
