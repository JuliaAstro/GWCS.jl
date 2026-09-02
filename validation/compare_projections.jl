using GWCS
using JSON
using StaticArrays

const PROJECTIONS = Dict(
    "TAN" => GnomonicProjection(),
    "STG" => StereographicProjection(),
    "ARC" => ZenithalEquidistant(),
    "ZEA" => ZenithalEqualArea(),
    "SIN" => SlantOrthographic(0.1, 0.2),
    "AZP" => ZenithalPerspective(1.5, 15.0),
    "SZP" => SlantZenithalPerspective(1.5, 15.0, 80.0),
    "AIR" => AiryProjection(75.8)
)

function main()
    python = get(ENV, "GWCS_PYTHON", joinpath(homedir(), ".venvs", "311_main", "bin", "python"))
    script = joinpath(@__DIR__, "python_projection_reference.py")
    records = JSON.parse(read(`$python $script`, String))
    max_forward_error = 0.0
    max_inverse_error = 0.0

    for record in records
        projection = PROJECTIONS[record["name"]]
        to_native = IntermediateToNative(projection)
        to_intermediate = NativeToIntermediate(projection)
        for sample in record["intermediate"]
            actual = to_native(SVector{2,Float64}(sample["input"]))
            expected = SVector{2,Float64}(sample["output"])
            longitude_error = abs(mod(actual[1] - expected[1] + 180, 360) - 180)
            abs(abs(expected[2]) - 90) <= 1e-12 && (longitude_error = 0.0)
            max_forward_error = max(max_forward_error, longitude_error, abs(actual[2] - expected[2]))
        end
        for sample in record["native"]
            actual = to_intermediate(SVector{2,Float64}(sample["input"]))
            expected = SVector{2,Float64}(sample["output"])
            max_inverse_error = max(max_inverse_error, maximum(abs, actual - expected))
        end
    end

    println("projections compared: ", length(records))
    println("maximum intermediate-to-native error (degrees): ", max_forward_error)
    println("maximum native-to-intermediate error (degrees): ", max_inverse_error)

    max_forward_error <= 1e-10 || error("intermediate-to-native comparison exceeded 1e-10 degrees")
    max_inverse_error <= 1e-10 || error("native-to-intermediate comparison exceeded 1e-10 degrees")
end

main()
