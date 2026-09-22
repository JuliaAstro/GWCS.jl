using GWCS
using JSON
using StaticArrays

include("python.jl")

function main()
    python = validation_python()
    reference_script = joinpath(@__DIR__, "python_bounding_box_reference.py")
    reference = JSON.parse(read(`$python $reference_script`, String))

    input = CoordinateFrame{3}(name="input", axes_names=("x", "y", "slice"))
    output = CoordinateFrame{3}(name="output")
    simple_box = BoundingBox{3}(((0.0, 5.0), (10.0, 15.0)), (1, 2))
    simple = WCS(Step(input, IdentityTransform{3}()), Step(output); bounding_box=simple_box)
    boxes = Dict((1,) => BoundingBox{3}(((0.0, 5.0), (0.0, 5.0)), (1, 2)),
                 (2,) => BoundingBox{3}(((10.0, 15.0), (10.0, 15.0)), (1, 2)))
    compound = WCS(Step(input, IdentityTransform{3}()), Step(output); bounding_box=CompoundBoundingBox(boxes, (3,)))
    scaled = parallel(Scale(2.0), Scale(3.0), Scale(1.0))
    dense = WCS(Step(input, scaled), Step(output); bounding_box=BoundingBox(((0.0, 5.0), (10.0, 15.0), (0.0, 200.0))))
    compound_mask = WCS(Step(input, scaled), Step(output); bounding_box=compound.bounding_box)

    @assert pixel_to_world(simple, 2.0, 12.0, 100.0) == SVector{3,Float64}(reference["simple_inside"])
    @assert all(isnan, pixel_to_world(simple, 7.0, 12.0, 100.0))
    @assert all(isnothing, reference["simple_outside"])
    @assert pixel_to_world(simple, 7.0, 12.0, 100.0; fill_value=-1.0) == SVector{3,Float64}(reference["simple_fill"])
    @assert pixel_to_world(compound, 2.0, 3.0, 1.0) == SVector{3,Float64}(reference["compound_inside"])
    @assert all(isnan, pixel_to_world(compound, 2.0, 3.0, 2.0))
    @assert all(isnothing, reference["compound_outside"])
    @assert all(isnan, pixel_to_world(compound, 2.0, 3.0, 3.0))
    @assert reference["missing_selector"] == "RuntimeError"
    @assert [in_image(dense, 4.0, 36.0, 100.0),
             in_image(dense, 14.0, 36.0, 100.0),
             in_image(dense, 14.0, 36.0, 100.0; with_bounding_box=false)] == reference["dense_in_image"]
    simple_worlds = [4.0 14.0 4.0; 36.0 36.0 36.0; 100.0 100.0 NaN]
    @assert collect(in_image(dense, simple_worlds)) == reference["dense_in_image_batch"]
    @assert !in_image(simple, 7.0, 12.0, 100.0)
    @assert reference["ignored_in_image_outside"]
    @assert [in_image(compound_mask, 4.0, 9.0, 1.0), in_image(compound_mask, 4.0, 9.0, 2.0)] == [true, false]
    @assert reference["compound_in_image"] == ["AttributeError", "AttributeError"]
    @assert !in_image(compound_mask, 4.0, 9.0, 3.0)
    @assert reference["compound_missing_in_image"] == "AttributeError"

    simple_grid = grid_from_bounding_box(BoundingBox(((-1.0, 2.9), (6.0, 7.5))); step=(1.0, 0.5), center=false)
    compound_grid = grid_from_bounding_box(compound.bounding_box; selector=1)
    @assert simple_grid == reduce(vcat, permutedims.(reference["simple_grid"]))
    @assert compound_grid[1:2, :] == reduce(vcat, permutedims.(reference["compound_grid"]))

    println("bounding-box scalar and mask cases compared: 15")
    println("bounding-box grid points compared: ", size(simple_grid, 2) + size(compound_grid, 2))
    println("Python unknown-selector behavior: RuntimeError; GWCS.jl treats it as outside the domain")
    println("Python in_image does not enforce ignored-input boxes and errors on compound boxes; GWCS.jl handles both")
end

main()
