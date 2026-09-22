using GWCS
using JSON
using StaticArrays

include("python.jl")

function parser_transforms()
    metadata = (
        identity="!<tag:stsci.edu:asdf/transform/identity-1.4.0> {n_dims: 2}",
        affine="""
        !<tag:stsci.edu:asdf/transform/affine-1.5.0>
        matrix: [[1.0, 0.2], [-0.1, 2.0]]
        translation: [3.0, 4.0]
        """,
        gnomonic="!<tag:stsci.edu:asdf/transform/gnomonic-1.4.0> {direction: pix2sky}",
        native2celestial="!<tag:stsci.edu:asdf/transform/rotate3d-1.5.0> {phi: 10.0, theta: 20.0, psi: 180.0, direction: native2celestial}",
        celestial2native="!<tag:stsci.edu:asdf/transform/rotate3d-1.5.0> {phi: 10.0, theta: 20.0, psi: 180.0, direction: celestial2native}",
        euler="!<tag:stsci.edu:asdf/transform/rotate3d-1.5.0> {phi: 1.0, theta: 2.0, psi: 3.0, direction: xyz}",
        spherical_sequence="""
        !<tag:stsci.edu:asdf/transform/rotate_sequence_3d-1.2.0>
        angles: [1.0, 2.0, 3.0]
        axes_order: zyx
        rotation_type: spherical
        """
    )
    return Dict(String(name) => parse_transform(GWCS._load_tagged_yaml(value)) for (name, value) in pairs(metadata))
end

function main()
    python = validation_python()
    reference_script = joinpath(@__DIR__, "python_parser_reference.py")
    records = JSON.parse(read(`$python $reference_script`, String))
    transforms = parser_transforms()
    max_forward_error = 0.0
    max_inverse_error = 0.0

    for record in records
        transform = transforms[record["transform"]]
        input = SVector{2,Float64}(record["input"])
        python_output = SVector{2,Float64}(record["output"])
        python_roundtrip = SVector{2,Float64}(record["roundtrip"])
        julia_output = transform(input)
        julia_roundtrip = inverse_transform(transform)(python_output)
        max_forward_error = max(max_forward_error, maximum(abs, julia_output - python_output))
        max_inverse_error = max(max_inverse_error, maximum(abs, julia_roundtrip - python_roundtrip))
    end

    println("parser points compared: ", length(records))
    println("maximum forward error: ", max_forward_error)
    println("maximum inverse error: ", max_inverse_error)
    max_forward_error <= 1e-12 || error("forward parser comparison exceeded 1e-12")
    max_inverse_error <= 1e-12 || error("inverse parser comparison exceeded 1e-12")
end

main()
