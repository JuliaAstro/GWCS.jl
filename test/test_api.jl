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

@testset "Coordinates without a result and strict evaluation" begin
    sky_header(code) = Dict("NAXIS" => 2, "CTYPE1" => "RA---$code", "CTYPE2" => "DEC--$code", "CRPIX1" => 181.0, "CRPIX2" => 91.0,
        "CRVAL1" => 0.0, "CRVAL2" => 0.0, "CDELT1" => -1.0, "CDELT2" => 1.0)
    ait = parse_fits_wcs(sky_header("AIT"))
    grid = [1.0 181.0 360.0; 1.0 91.0 180.0]
    # By default, off-sky pixels of an all-sky projection return fill_value.
    filled = pixel_to_world(ait, grid)
    @test all(isnan, filled[:, 1]) && all(isnan, filled[:, 3])
    @test filled[:, 2] ≈ [0.0, 0.0] atol=1e-12
    @test pixel_to_world(ait, grid; fill_value=-999.0)[:, 1] == [-999.0, -999.0]
    @test all(isnan, pixel_to_world(ait, 1.0, 1.0))
    @test_throws DomainError pixel_to_world(ait, grid; strict=true)
    @test @inferred(pixel_to_world(ait, SVector(181.0, 91.0); strict=true)) == pixel_to_world(ait, SVector(181.0, 91.0))

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
    @test_throws r"pixel \(1\.0, 1\.0\) has no world coordinate: IntermediateToNative\(HammerAitoffProjection\) in the step from frame \"intermediate\" to \"world\" is undefined at" pixel_to_world(ait, 1.0, 1.0; strict=true)
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
