@testset "Bounding-box domains" begin
    dense = BoundingBox(((-0.5, 9.5), (-0.5, 19.5)))
    @test input_dimension(dense) == 2
    @test bounded_axes(dense) == (1, 2)
    @test ignored_axes(dense) == ()
    @test bounds(dense, 1) == (-0.5, 9.5)
    @test dense isa BoundingBox{2,Float64}
    @test BoundingBox{2,Float64}(((-0.5, 9.5), (-0.5, 19.5)), :F).intervals == dense.intervals
    @test SVector(1.0, 2.0) in dense
    @test !(SVector(10.0, 2.0) in dense)

    c_order = BoundingBox(((10.0, 20.0), (1.0, 2.0)); order=:C)
    @test c_order.intervals == ((1.0, 2.0), (10.0, 20.0))
    @test SVector(1.5, 15.0) in c_order

    partial = BoundingBox{3}(((0.0, 5.0), (10.0, 20.0)), (1, 3))
    @test bounded_axes(partial) == (1, 3)
    @test ignored_axes(partial) == (2,)
    @test bounds(partial, 2) === nothing
    @test SVector(2.0, -100.0, 15.0) in partial
    @test !(SVector(6.0, -100.0, 15.0) in partial)

    @test_throws "exactly two bounds" BoundingBox(((0.0, 1.0, 2.0),))
    @test_throws "lower limits must not exceed" BoundingBox(((2.0, 1.0),))
    @test_throws "bounded axes must be unique" BoundingBox{2}(((0.0, 1.0), (2.0, 3.0)), (1, 1))

    boxes = Dict((1,) => BoundingBox{3}(((0.0, 5.0), (0.0, 5.0)), (1, 2)),
                 (2,) => BoundingBox{3}(((10.0, 15.0), (10.0, 15.0)), (1, 2)))
    compound = CompoundBoundingBox(boxes, (3,))
    @test selector_axes(compound) == (3,)
    @test selector_ignored(compound) == (true,)
    @test select_bounding_box(compound, 1).intervals == ((0.0, 5.0), (0.0, 5.0))
    @test compound[2].intervals == ((10.0, 15.0), (10.0, 15.0))
    @test SVector(2.0, 3.0, 1.0) in compound
    @test !(SVector(7.0, 3.0, 1.0) in compound)
    @test !(SVector(2.0, 3.0, 3.0) in compound)
    @test_throws KeyError compound[3]
    @test_throws DimensionMismatch select_bounding_box(CompoundBoundingBox(Dict((1, 2) => BoundingBox{3}(((0.0, 1.0),), (3,))), (1, 2)), 1)
    @test_throws "at least one entry" CompoundBoundingBox(Dict(), (1,))
    @test_throws "flags must be booleans" CompoundBoundingBox(boxes, (3,); selector_ignored=(1,))

    bounded_selector = CompoundBoundingBox(Dict((2,) => BoundingBox{2}(((1.0, 3.0), (10.0, 12.0)), (1, 2))), (1,); selector_ignored=false)
    @test SVector(2.0, 11.0) in bounded_selector
    @test !(SVector(4.0, 11.0) in bounded_selector)
    @test_throws "ignore flag is inconsistent" CompoundBoundingBox(boxes, (3,); selector_ignored=false)

    grid = grid_from_bounding_box(BoundingBox(((-1.0, 2.9), (6.0, 7.5))); step=(1.0, 0.5), center=false)
    @test size(grid) == (2, 20)
    @test grid[:, 1] == [-1.0, 6.0]
    @test grid[:, end] == [3.0, 7.5]

    centered = grid_from_bounding_box(BoundingBox(((-1.0, 2.9),)))
    @test centered == reshape([-1.0, 0.0, 1.0, 2.0, 3.0], 1, :)

    partial_grid = grid_from_bounding_box(partial; step=(1.0, 5.0), fixed=(2 => 7.0,))
    @test size(partial_grid) == (3, 18)
    @test all(==(7.0), partial_grid[2, :])
    @test_throws "fixed values are required for ignored axes" grid_from_bounding_box(partial)
    @test_throws "grid step has length" grid_from_bounding_box(dense; step=(1.0,))

    compound_grid = grid_from_bounding_box(compound; selector=1)
    @test size(compound_grid) == (3, 36)
    @test all(==(1.0), compound_grid[3, :])
    @test_throws "selector is required" grid_from_bounding_box(compound)
end

@testset "WCS compound domains" begin
    input = CoordinateFrame{3}(name="input", axes_names=("x", "y", "slice"))
    output = CoordinateFrame{3}(name="output")
    boxes = Dict((1,) => BoundingBox{3}(((0.0, 5.0), (0.0, 5.0)), (1, 2)),
                 (2,) => BoundingBox{3}(((10.0, 15.0), (10.0, 15.0)), (1, 2)))
    compound = CompoundBoundingBox(boxes, (3,))
    wcs = WCS(Step(input, IdentityTransform{3}()), Step(output); bounding_box=compound, pixel_shape=(20, 20, 2))

    @test pixel_to_world(wcs, 2.0, 3.0, 1.0) == SVector(2.0, 3.0, 1.0)
    @test all(isnan, pixel_to_world(wcs, 7.0, 3.0, 1.0))
    @test pixel_to_world(wcs, 7.0, 3.0, 1.0; with_bounding_box=false) == SVector(7.0, 3.0, 1.0)
    @test pixel_to_world(wcs, 7.0, 3.0, 1.0; fill_value=-1.0) == SVector(-1.0, -1.0, -1.0)
    @test all(isnan, pixel_to_world(wcs, 2.0, 3.0, 3.0))
    @test world_to_pixel(wcs, 2.0, 3.0, 1.0) == SVector(2.0, 3.0, 1.0)
    @test all(isnan, world_to_pixel(wcs, 7.0, 3.0, 1.0))

    pixels = [2.0 7.0 12.0; 3.0 3.0 12.0; 1.0 1.0 2.0]
    worlds = pixel_to_world(wcs, pixels)
    @test worlds[:, 1] == pixels[:, 1]
    @test all(isnan, worlds[:, 2])
    @test worlds[:, 3] == pixels[:, 3]

    sliced = fix_inputs(wcs, "slice" => 1.0)
    @test sliced.bounding_box isa BoundingBox{2}
    @test sliced.bounding_box.intervals == ((0.0, 5.0), (0.0, 5.0))
    @test pixel_to_world(sliced, 2.0, 3.0) == SVector(2.0, 3.0, 1.0)
    @test_throws DomainError fix_inputs(wcs, "slice" => 3.0)

    bounded_slice = fix_inputs(wcs, "x" => 2.0)
    @test bounded_slice.bounding_box isa CompoundBoundingBox{2}
    @test Set(keys(bounded_slice.bounding_box.boxes)) == Set(((1,),))
    @test pixel_to_world(bounded_slice, 3.0, 1.0) == SVector(2.0, 3.0, 1.0)
    @test_throws DomainError fix_inputs(wcs, "x" => 7.0)

    replacement = BoundingBox(((0.0, 1.0), (0.0, 1.0), (1.0, 2.0)))
    replaced = with_bounding_box(wcs, replacement)
    @test replaced.bounding_box === replacement
    @test with_bounding_box(replaced, nothing).bounding_box === nothing
    @test with_pixel_shape(wcs, (30, 40, 2)).pixel_shape == (30, 40, 2)
    @test with_pixel_shape(wcs, nothing).pixel_shape === nothing

    grid = grid_from_bounding_box(wcs; selector=1)
    @test size(grid) == (3, 36)
    @test all(==(1.0), grid[3, :])

    allocated_domain(wcs, coords) = @allocated pixel_to_world(wcs, coords)
    coords = SVector(2.0, 3.0, 1.0)
    pixel_to_world(wcs, coords)
    @test @inferred(pixel_to_world(wcs, coords)) == coords
    @test allocated_domain(wcs, coords) == 0

    no_inverse = CallableTransform{3,3}((x, y, slice) -> (x, y, slice))
    numerical = WCS(Step(input, no_inverse), Step(output); bounding_box=compound)
    @test_throws "requires an initial_guess" world_to_pixel(numerical, 2.0, 3.0, 1.0)
    @test world_to_pixel(numerical, 2.0, 3.0, 1.0; initial_guess=(2.0, 3.0, 1.0)) ≈ SVector(2.0, 3.0, 1.0)

    input4 = CoordinateFrame{4}(name="input4", axes_names=("x", "y", "detector", "slice"))
    output4 = CoordinateFrame{4}(name="output4")
    entries = Dict((detector, slice) => BoundingBox{4}(((0.0, 5.0), (0.0, 5.0)), (1, 2)) for detector in (1, 2) for slice in (10, 20))
    compound4 = CompoundBoundingBox(entries, (3, 4))
    wcs4 = WCS(Step(input4, IdentityTransform{4}()), Step(output4); bounding_box=compound4)
    reduced = fix_inputs(wcs4, "detector" => 1)
    @test reduced.bounding_box isa CompoundBoundingBox{3}
    @test selector_axes(reduced.bounding_box) == (3,)
    @test Set(keys(reduced.bounding_box.boxes)) == Set(((10,), (20,)))
end
