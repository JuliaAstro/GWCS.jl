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

    batch = [2.0 3.0; 4.0 5.0]
    world_batch = pixel_to_world(wcs, batch)
    @test world_batch == [5.0 7.0; 11.0 14.0]
    @test world_to_pixel(wcs, world_batch) ≈ batch
    @test_throws DimensionMismatch pixel_to_world(wcs, [1.0])
    @test_throws DimensionMismatch world_to_pixel(wcs, [1.0])

    coords = SVector(2.0, 4.0)
    allocated_pixel_to_world(wcs, coords) = @allocated pixel_to_world(wcs, coords)
    allocated_world_to_pixel(wcs, coords) = @allocated world_to_pixel(wcs, coords)
    pixel_to_world(wcs, coords)
    world_to_pixel(wcs, SVector(5.0, 11.0))
    @test @inferred(pixel_to_world(wcs, coords)) == SVector(5.0, 11.0)
    @test @inferred(world_to_pixel(wcs, SVector(5.0, 11.0))) ≈ coords
    @test allocated_pixel_to_world(wcs, coords) == 0
    @test allocated_world_to_pixel(wcs, SVector(5.0, 11.0)) == 0
end
