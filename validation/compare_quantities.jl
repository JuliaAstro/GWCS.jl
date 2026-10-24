using ASDF
using GWCS
using JSON
using Logging

include("python.jl")

# Parameters read as ASDF quantities must be written back in their units, and Python must read GWCS.jl's rewrite of
# test/fixtures/quantities.asdf with the same frames, units, and results.
function main()
    fixture_path = joinpath(@__DIR__, "..", "test", "fixtures", "quantities.asdf")
    fixture = with_logger(NullLogger()) do
        ASDF.load(fixture_path; extensions=true)
    end
    wcss = Dict(name => parse_wcs(node) for (name, node) in fixture["wcs"])
    transforms = Dict(name => parse_transform(node) for (name, node) in fixture["transforms"])
    max_error = 0.0
    mktempdir() do directory
        filename = joinpath(directory, "quantities.asdf")
        write_asdf(filename, Dict("wcs" => wcss, "transforms" => transforms))
        records = JSON.parse(read(`$(validation_python()) $(joinpath(@__DIR__, "python_quantity_reference.py")) $filename $fixture_path`, String))
        for record in sort(records; by=record -> record["name"])
            max_error = max(max_error, record["error"])
            frame = haskey(record, "reference_frame") ? "  frame $(record["reference_frame"])$(record["frame_data"] ? " with data" : "")" : ""
            println(rpad(record["name"], 22), " error ", record["error"], frame, "  units ", join(("$name: $unit" for (name, unit) in sort(collect(record["units"]))), ", "))
        end
    end
    println("maximum error: ", max_error)
    max_error <= 1e-10 || error("quantity round-trip comparison exceeded 1e-10")
end

main()
