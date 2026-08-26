using GWCS
using JSON
using StaticArrays

function main()
    python = get(ENV, "GWCS_PYTHON", joinpath(homedir(), ".venvs", "311_main", "bin", "python"))
    reference_script = joinpath(@__DIR__, "python_bounding_box_reference.py")
    reference = JSON.parse(read(`$python $reference_script`, String))

    input = CoordinateFrame{3}(name="input", axes_names=("x", "y", "slice"))
    output = CoordinateFrame{3}(name="output")
    simple_box = BoundingBox{3}(((0.0, 5.0), (10.0, 15.0)), (1, 2))
    simple = WCS(Step(input, IdentityTransform{3}()), Step(output); bounding_box=simple_box)
    boxes = Dict((1,) => BoundingBox{3}(((0.0, 5.0), (0.0, 5.0)), (1, 2)),
                 (2,) => BoundingBox{3}(((10.0, 15.0), (10.0, 15.0)), (1, 2)))
    compound = WCS(Step(input, IdentityTransform{3}()), Step(output); bounding_box=CompoundBoundingBox(boxes, (3,)))

    @assert pixel_to_world(simple, 2.0, 12.0, 100.0) == SVector{3,Float64}(reference["simple_inside"])
    @assert all(isnan, pixel_to_world(simple, 7.0, 12.0, 100.0))
    @assert all(isnothing, reference["simple_outside"])
    @assert pixel_to_world(simple, 7.0, 12.0, 100.0; fill_value=-1.0) == SVector{3,Float64}(reference["simple_fill"])
    @assert pixel_to_world(compound, 2.0, 3.0, 1.0) == SVector{3,Float64}(reference["compound_inside"])
    @assert all(isnan, pixel_to_world(compound, 2.0, 3.0, 2.0))
    @assert all(isnothing, reference["compound_outside"])
    @assert all(isnan, pixel_to_world(compound, 2.0, 3.0, 3.0))
    @assert reference["missing_selector"] == "RuntimeError"

    simple_grid = grid_from_bounding_box(BoundingBox(((-1.0, 2.9), (6.0, 7.5))); step=(1.0, 0.5), center=false)
    compound_grid = grid_from_bounding_box(compound.bounding_box; selector=1)
    @assert simple_grid == reduce(vcat, permutedims.(reference["simple_grid"]))
    @assert compound_grid[1:2, :] == reduce(vcat, permutedims.(reference["compound_grid"]))

    println("bounding-box scalar cases compared: 6")
    println("bounding-box grid points compared: ", size(simple_grid, 2) + size(compound_grid, 2))
    println("Python unknown-selector behavior: RuntimeError; GWCS.jl treats it as outside the domain")
end

main()
