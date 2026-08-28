module RegistrationFixture

using GWCS
using StaticArrays
import GWCS: compile_transform, inverse_transform

const parsed_context = Ref{Any}()
const validated_context = Ref{Any}()

struct OffsetTransform{T<:Real} <: AbstractTransform{1,1}
    offset::T
end

@inline function (transform::OffsetTransform)(coords::StaticVector{1})
    return SVector(coords[1] + transform.offset)
end

inverse_transform(transform::OffsetTransform) = OffsetTransform(-transform.offset)
compile_transform(transform::OffsetTransform) = Shift(transform.offset)

struct WrappedTransform{T<:AbstractTransform{1,1}} <: AbstractTransform{1,1}
    transform::T
end

@inline (transform::WrappedTransform)(coords::StaticVector{1}) = transform.transform(coords)
inverse_transform(transform::WrappedTransform) = WrappedTransform(inverse_transform(transform.transform))
compile_transform(transform::WrappedTransform) = compile_transform(transform.transform)

struct BadCompiled <: AbstractTransform{1,1} end
@inline (::BadCompiled)(coords::StaticVector{1}) = coords
compile_transform(::BadCompiled) = IdentityTransform{2}()

struct BadInverse <: AbstractTransform{1,1} end
@inline (::BadInverse)(coords::StaticVector{1}) = coords
inverse_transform(::BadInverse) = IdentityTransform{2}()

function parse_offset(node, context)
    parsed_context[] = (metadata_tag(context), metadata_version(context), metadata_path(context))
    offset = metadata_value(required_metadata(node, "offset", context))
    return OffsetTransform(offset)
end

function parse_wrapped(node, context)
    child = required_metadata(node, "transform", context)
    transform = parse_transform(child, context; path="$(metadata_path(context))/transform")
    return WrappedTransform(transform)
end

parse_default(node, context) = OffsetTransform(metadata_value(required_metadata(node, "offset", context)))
parse_other(node, context) = OffsetTransform(10 * metadata_value(required_metadata(node, "offset", context)))

function validate_legacy(node, context)
    validated_context[] = (metadata_tag(context), metadata_version(context), metadata_path(context))
    required_metadata(node, "legacy_offset", context)
    return nothing
end

function normalize_legacy(node, context)
    normalized = Dict{Any,Any}(node)
    normalized["offset"] = metadata_value(node["legacy_offset"])
    return normalized
end

end

@testset "Transform registration" begin
    registry = copy(default_transform_registry())
    register_transform!(registry, "tag:example.org:gwcs/offset", RegistrationFixture.parse_offset; min_version=v"1.0.0", max_version=v"2.0.0")
    register_transform!(registry, "tag:example.org:gwcs/wrapped", RegistrationFixture.parse_wrapped; min_version=v"1.0.0", max_version=v"2.0.0")
    register_transform!(registry, "tag:other.example:gwcs/offset", RegistrationFixture.parse_other; min_version=v"1.0.0", max_version=v"2.0.0")

    offset_metadata = """
    !<tag:example.org:gwcs/offset-1.2.3>
    offset: 2.5
    """
    transform = parse_transform(GWCS._load_tagged_yaml(offset_metadata); registry)
    @test transform isa RegistrationFixture.OffsetTransform
    @test input_dimension(transform) == 1
    @test output_dimension(transform) == 1
    @test transform(SVector(1.0)) == SVector(3.5)
    @test inverse_transform(transform)(SVector(3.5)) == SVector(1.0)
    @test compile_transform(transform) == Shift(2.5)
    @test RegistrationFixture.parsed_context[] == ("tag:example.org:gwcs/offset", v"1.2.3", "transform")
    @test materialize_metadata([1, 2]) == [1, 2]

    annotated_metadata = """
    !<tag:example.org:gwcs/offset-1.2.3>
    offset: 2.5
    name: registered
    inputs: [x]
    outputs: [u]
    """
    annotated = parse_transform(GWCS._load_tagged_yaml(annotated_metadata); registry)
    @test annotated isa AnnotatedTransform
    @test annotated.transform isa RegistrationFixture.OffsetTransform
    @test transform_name(annotated) == "registered"
    @test input_names(annotated) == ("x",)
    @test output_names(annotated) == ("u",)

    register_transform!(registry, "tag:example.org:gwcs/preannotated", (node, context) -> AnnotatedTransform(RegistrationFixture.OffsetTransform(1.0); inputs=("x",), outputs=("u",)))
    partially_annotated = parse_transform(GWCS._load_tagged_yaml("""
    !<tag:example.org:gwcs/preannotated-1.0.0>
    name: retained_axes
    """); registry)
    @test transform_name(partially_annotated) == "retained_axes"
    @test input_names(partially_annotated) == ("x",)
    @test output_names(partially_annotated) == ("u",)

    wrapped_metadata = """
    !<tag:example.org:gwcs/wrapped-1.0.0>
    transform: !<tag:example.org:gwcs/offset-1.9.9>
      offset: 4.0
    """
    wrapped = parse_transform(GWCS._load_tagged_yaml(wrapped_metadata); registry)
    @test wrapped isa RegistrationFixture.WrappedTransform
    @test wrapped(SVector(1.0)) == SVector(5.0)

    other_metadata = replace(offset_metadata, "tag:example.org" => "tag:other.example")
    @test parse_transform(GWCS._load_tagged_yaml(other_metadata); registry)(SVector(1.0)) == SVector(26.0)

    explicit_inverse_metadata = """
    !<tag:example.org:gwcs/offset-1.0.0>
    offset: 2.0
    inverse: !<tag:example.org:gwcs/offset-1.0.0>
      offset: -2.0
    """
    explicit = parse_transform(GWCS._load_tagged_yaml(explicit_inverse_metadata); registry)
    @test explicit isa WithInverse
    @test inverse_transform(explicit)(explicit(SVector(3.0))) == SVector(3.0)

    annotated_inverse_metadata = """
    !<tag:example.org:gwcs/offset-1.0.0>
    offset: 2.0
    name: named_inverse
    inputs: [x]
    outputs: [u]
    inverse: !<tag:example.org:gwcs/offset-1.0.0>
      offset: -2.0
    """
    annotated_explicit = parse_transform(GWCS._load_tagged_yaml(annotated_inverse_metadata); registry)
    annotated_backward = inverse_transform(annotated_explicit)
    @test input_names(annotated_backward) == ("u",)
    @test output_names(annotated_backward) == ("x",)
    @test annotated_backward(annotated_explicit(SVector(3.0))) == SVector(3.0)

    unsupported_metadata = replace(offset_metadata, "offset-1.2.3" => "offset-2.0.0")
    @test_throws "registered ranges are [1.0.0, 2.0.0)" parse_transform(GWCS._load_tagged_yaml(unsupported_metadata); registry)
    @test_throws "overlaps existing range" register_transform!(registry, "tag:example.org:gwcs/offset", RegistrationFixture.parse_offset; min_version=v"1.5.0", max_version=v"2.5.0")
    @test_throws "must omit the schema version" register_transform!(registry, "tag:example.org:gwcs/versioned-1.0.0", RegistrationFixture.parse_offset)

    versioned_registry = copy(default_transform_registry())
    register_transform!(versioned_registry, "tag:example.org:gwcs/versioned", RegistrationFixture.parse_offset; min_version=v"1.0.0", max_version=v"2.0.0")
    register_transform!(versioned_registry, "tag:example.org:gwcs/versioned", RegistrationFixture.parse_other; min_version=v"2.0.0", max_version=v"3.0.0")
    versioned_metadata = "!<tag:example.org:gwcs/versioned-2.0.0> {offset: 2.0}"
    @test parse_transform(GWCS._load_tagged_yaml(versioned_metadata); registry=versioned_registry)(SVector(1.0)) == SVector(21.0)
    unsupported_version = replace(versioned_metadata, "versioned-2.0.0" => "versioned-3.0.0")
    @test_throws "unsupported transform tag version" parse_transform(GWCS._load_tagged_yaml(unsupported_version); registry=versioned_registry)

    exact_registry = copy(default_transform_registry())
    register_transform!(exact_registry, "tag:example.org:gwcs/exact", RegistrationFixture.parse_offset; versions=(v"1.0.0", v"1.2.0"), validate=RegistrationFixture.validate_legacy, normalize=RegistrationFixture.normalize_legacy)
    exact_metadata = "!<tag:example.org:gwcs/exact-1.2.0> {legacy_offset: 2.0}"
    @test parse_transform(GWCS._load_tagged_yaml(exact_metadata); registry=exact_registry)(SVector(1.0)) == SVector(3.0)
    @test RegistrationFixture.validated_context[] == ("tag:example.org:gwcs/exact", v"1.2.0", "transform")
    @test_throws "supported versions are {1.0.0, 1.2.0}" parse_transform(GWCS._load_tagged_yaml(replace(exact_metadata, "exact-1.2.0" => "exact-1.1.0")); registry=exact_registry)
    @test_throws "cannot be combined" register_transform!(exact_registry, "tag:example.org:gwcs/mixed", RegistrationFixture.parse_offset; versions=(v"1.0.0",), min_version=v"1.0.0")

    fresh_registry = TransformRegistry()
    register_transform!(fresh_registry, "tag:example.org:gwcs/offset", RegistrationFixture.parse_offset; versions=(v"1.0.0",))
    fresh_wcs = parse_wcs("""
    !<tag:stsci.edu:gwcs/wcs-1.4.0>
    name: custom
    steps:
      - !<tag:stsci.edu:gwcs/step-1.3.0>
        frame: !<tag:stsci.edu:gwcs/frame-1.2.0> {name: input, naxes: 1}
        transform: !<tag:example.org:gwcs/offset-1.0.0> {offset: 2.0}
      - !<tag:stsci.edu:gwcs/step-1.3.0>
        frame: !<tag:stsci.edu:gwcs/frame-1.2.0> {name: output, naxes: 1}
    """; registry=fresh_registry)
    @test pixel_to_world(fresh_wcs, 1.0) == SVector(3.0)

    bad_registry = copy(default_transform_registry())
    register_transform!(bad_registry, "tag:example.org:gwcs/bad", (node, context) -> 1)
    bad_metadata = "!<tag:example.org:gwcs/bad-1.0.0> {}"
    @test_throws "expected AbstractTransform" parse_transform(GWCS._load_tagged_yaml(bad_metadata); registry=bad_registry)

    frame = CoordinateFrame{1}(; name="input")
    output = CoordinateFrame{1}(; name="output")
    @test_throws "compiled form" WCS(Step(frame, RegistrationFixture.BadCompiled()), Step(output))
    @test_throws "inverse of" WCS(Step(frame, RegistrationFixture.BadInverse()), Step(output))

    default_tag = "tag:example.org:gwcs/default_offset_$(objectid(registry))"
    register_transform!(default_tag, RegistrationFixture.parse_default; min_version=v"1.0.0", max_version=v"1.1.0")
    default_metadata = "!<$default_tag-1.0.0> {offset: 3.0}"
    @test parse_transform(GWCS._load_tagged_yaml(default_metadata))(SVector(2.0)) == SVector(5.0)

    roman_metadata = read(joinpath(@__DIR__, "fixtures", "roman_wcs.yaml"), String)
    custom_roman_metadata = replace(roman_metadata, "!transform/shift-1.4.0 {offset: -10.0}" => "!<tag:example.org:gwcs/offset-1.2.3> {offset: -10.0}"; count=1)
    custom_wcs = parse_wcs(custom_roman_metadata; registry)
    @test pixel_to_world(custom_wcs, 11.0, 22.0) ≈ SVector(1.0, 2.0) atol=1e-12
    @test world_to_pixel(custom_wcs, 1.0, 2.0) ≈ SVector(11.0, 22.0) atol=1e-10

    compiled_wcs = WCS(Step(frame, transform), Step(output))
    @test compiled_wcs.forward isa Shift
    @test pixel_to_world(compiled_wcs, 1.0) == SVector(3.5)
    @test world_to_pixel(compiled_wcs, 3.5) == SVector(1.0)
    allocated_forward(wcs, coords) = @allocated pixel_to_world(wcs, coords)
    allocated_inverse(wcs, coords) = @allocated world_to_pixel(wcs, coords)
    pixel_to_world(compiled_wcs, SVector(1.0))
    world_to_pixel(compiled_wcs, SVector(3.5))
    @test @inferred(pixel_to_world(compiled_wcs, SVector(1.0))) == SVector(3.5)
    @test @inferred(world_to_pixel(compiled_wcs, SVector(3.5))) == SVector(1.0)
    @test allocated_forward(compiled_wcs, SVector(1.0)) == 0
    @test allocated_inverse(compiled_wcs, SVector(3.5)) == 0
end
