using ASDF
using GWCS
using GWCS: AffineTransform, ArrayLabelMapper, Polynomial1D, Polynomial2D, RegionSelector, RemapAxes, Rotation2D, Scale, Shift, compose,
            input_dimension, parallel, transform
using JSON
using StaticArrays

include("python.jl")

# Canonical forms that compose() and parallel() build, and every array GWCS.jl writes as an ndarray, must stay readable by
# Python GWCS.
function main()
    distortion = (Polynomial2D([0.1 0.2 0.3; 1.0 0.01 0.0; 0.02 0.0 0.0]; domain=((0.0, 10.0), (-5.0, 5.0)), window=((-1.0, 1.0), (-1.0, 1.0))), Polynomial2D([0.5 1.1 0.03; -0.2 0.0 0.0; 0.04 0.0 0.0]))
    linear = (Polynomial2D([0.1 0.2; 1.0 0.0]; domain=((0.0, 10.0), (-5.0, 5.0)), window=((-1.0, 1.0), (-1.0, 1.0))), Polynomial2D([0.5 1.1; -0.2 0.0]))
    transforms = Dict(
        "three_shifts" => Shift(1.0) & Shift(2.0) & Shift(3.0),
        "rotation_chain" => Shift(1.0) & Shift(2.0) | Rotation2D(-75.0) | Scale(2.0) & Scale(0.5),
        "affine_blocks" => AffineTransform([1.0 0.2; -0.1 2.0], [3.0, 4.0]) & Scale(2.0),
        "polynomial_pair" => compose(RemapAxes{2}((1, 2, 1, 2)), parallel(distortion...)),
        "linear_pair" => compose(AffineTransform([2.0 0.0; 0.0 3.0], [1.0, 1.0]), RemapAxes{2}((1, 2, 1, 2)), parallel(linear...), Shift(1.0) & Shift(2.0)),
        "polynomial_1d" => Polynomial1D([0.5, 2.0, 0.01]) & Polynomial1D([-1.0, 0.5]),
        "region_selector" => RegionSelector(ArrayLabelMapper(Int16[0 2; 1 0]; no_label=0), (Int16(1) => Shift(1.0) & Shift(2.0), Int16(2) => Scale(2.0) & Scale(3.0)))
    )
    inputs = Dict(name => [collect(range(-2.0, 3.0; length=input_dimension(transform))) .+ offset for offset in (0.0, 0.5, 4.0)] for (name, transform) in transforms)
    # Positions inside the label map, where both readers select a branch.
    inputs["region_selector"] = [[0.0, 1.0], [1.0, 0.0], [0.2, 0.9]]
    max_error = 0.0
    mktempdir() do directory
        filename = joinpath(directory, "transforms.asdf")
        ASDF.write_file(filename, Dict("transforms" => transforms, "inputs" => inputs))
        records = JSON.parse(read(`$(validation_python()) $(joinpath(@__DIR__, "python_written_asdf_reference.py")) $filename`, String))
        length(records) == length(transforms) || error("Python read $(length(records)) of $(length(transforms)) transforms")
        for record in records
            transform = transforms[record["transform"]]
            for (coords, expected) in zip(inputs[record["transform"]], record["outputs"])
                max_error = max(max_error, maximum(abs, transform(SVector{length(coords)}(coords)) .- expected))
            end
            println(rpad(record["transform"], 16), " Python reads ", record["type"])
        end
    end
    println("maximum forward error: ", max_error)
    max_error <= 1e-12 || error("written ASDF comparison exceeded 1e-12")
end

main()
