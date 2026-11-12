using GWCS
using GWCS: parse_transform
using JSON
using StaticArrays

include("python.jl")

function main()
    python = validation_python()
    reference_script = joinpath(@__DIR__, "python_arithmetic_reference.py")
    records = JSON.parse(read(`$python $reference_script`, String))
    fixture = GWCS._load_tagged_yaml(read(joinpath(@__DIR__, "..", "test", "fixtures", "asdf_arithmetic.yaml"), String))
    transforms = Dict(name => parse_transform(node) for (name, node) in fixture["arithmetic"])
    max_error = 0.0

    for record in records
        transform = transforms[record["transform"]]
        result = only(transform(SVector(record["input"])))
        max_error = max(max_error, abs(result - record["output"]))
    end

    println("arithmetic points compared: ", length(records))
    println("maximum forward error: ", max_error)
    max_error <= 1e-12 || error("arithmetic comparison exceeded 1e-12")
end

main()
