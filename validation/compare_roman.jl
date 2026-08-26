using ASDF
using GWCS
using JSON
using Logging
using StaticArrays

function main(args)
    isempty(args) && error("usage: julia --project=validation validation/compare_roman.jl FILE...")

    python = get(ENV, "GWCS_PYTHON", joinpath(homedir(), ".venvs", "311_main", "bin", "python"))
    reference_script = joinpath(@__DIR__, "python_reference.py")
    records = JSON.parse(read(Cmd(vcat([python, reference_script], args)), String))
    wcs_cache = Dict{String,Any}()
    max_forward_error = 0.0
    max_inverse_error = 0.0

    for record in records
        filename = record["file"]
        wcs = get!(wcs_cache, filename) do
            asdf = with_logger(NullLogger()) do
                ASDF.load(filename; extensions=true, validate_checksum=false)
            end
            parse_wcs(asdf["roman"]["meta"]["wcs"])
        end

        pixel = SVector{2,Float64}(record["pixel"])
        python_world = SVector{2,Float64}(record["world"])
        python_roundtrip = SVector{2,Float64}(record["roundtrip"])
        julia_world = pixel_to_world(wcs, pixel)
        julia_inverse = world_to_pixel(wcs, python_world)
        max_forward_error = max(max_forward_error, maximum(abs, julia_world - python_world))
        max_inverse_error = max(max_inverse_error, maximum(abs, julia_inverse - python_roundtrip))
    end

    println("files compared: ", length(wcs_cache))
    println("points compared: ", length(records))
    println("maximum forward error (degrees): ", max_forward_error)
    println("maximum inverse error (pixels): ", max_inverse_error)

    max_forward_error <= 1e-10 || error("forward comparison exceeded 1e-10 degrees")
    max_inverse_error <= 1e-8 || error("inverse comparison exceeded 1e-8 pixels")
end

main(ARGS)
