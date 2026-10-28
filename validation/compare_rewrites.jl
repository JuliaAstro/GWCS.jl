using ASDF
using GWCS
using JSON
using Logging

include("python.jl")

# Python must read GWCS.jl's rewrite of each Python-written fixture with the same frames, axis metadata, parameter units,
# and results:
# quantities.asdf covers quantity parameters and astropy frames, projections.asdf every non-zenithal projection tag, and
# axes_order.asdf frames that reorder their axes, and metadata_fidelity.asdf parameter constraints and frame2d frames.
function compare(fixture_name)
    println(fixture_name)
    fixture_path = joinpath(@__DIR__, "..", "test", "fixtures", fixture_name)
    fixture = with_logger(NullLogger()) do
        ASDF.load(fixture_path; extensions=true)
    end
    wcss = Dict(name => parse_wcs(node) for (name, node) in fixture["wcs"])
    transforms = Dict(name => parse_transform(node) for (name, node) in fixture["transforms"])
    max_error = 0.0
    mktempdir() do directory
        filename = joinpath(directory, fixture_name)
        write_asdf(filename, Dict("wcs" => wcss, "transforms" => transforms))
        records = JSON.parse(read(`$(validation_python()) $(joinpath(@__DIR__, "python_rewrite_reference.py")) $filename $fixture_path`, String))
        for record in sort(records; by=record -> record["name"])
            max_error = max(max_error, record["error"])
            get(record, "axes_match", true) || error("Python reads $(record["name"]) of GWCS.jl's rewrite of $fixture_name with different axis metadata")
            get(record, "frames_match", true) || error("Python reads $(record["name"]) of GWCS.jl's rewrite of $fixture_name with other frame types")
            get(record, "constraints_match", true) || error("Python reads $(record["name"]) of GWCS.jl's rewrite of $fixture_name with other parameter constraints")
            frame = haskey(record, "reference_frame") ? "  frame $(record["reference_frame"])$(record["frame_data"] ? " with data" : "")" : ""
            println(rpad(record["name"], 22), " error ", record["error"], frame, "  units ", join(("$name: $unit" for (name, unit) in sort(collect(record["units"]))), ", "))
        end
    end
    println("maximum error: ", max_error)
    max_error <= 1e-10 || error("$fixture_name round-trip comparison exceeded 1e-10")
end

foreach(compare, ("quantities.asdf", "projections.asdf", "axes_order.asdf", "metadata_fidelity.asdf"))
