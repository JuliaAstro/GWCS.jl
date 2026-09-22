module AnnotationFixture

using GWCS
using StaticArrays
import GWCS: compile_transform

struct OffsetTransform{T<:Real} <: AbstractTransform{1,1}
    offset::T
end

@inline (transform::OffsetTransform)(coords::StaticVector{1}) = SVector(coords[1] + transform.offset)
compile_transform(transform::OffsetTransform) = Shift(transform.offset)

end

@testset "Transform annotations" begin
    annotated = AnnotatedTransform(Shift(2.0); name="offset", inputs=("x",), outputs=("u",))
    @test transform_name(annotated) == "offset"
    @test input_names(annotated) == ("x",)
    @test output_names(annotated) == ("u",)
    @test annotated(SVector(1.0)) == SVector(3.0)
    @test transform_jacobian(annotated, SVector(1.0)) == SMatrix{1,1}(1.0)

    inverse = inverse_transform(annotated)
    @test transform_name(inverse) == "offset"
    @test input_names(inverse) == ("u",)
    @test output_names(inverse) == ("x",)
    @test inverse(annotated(SVector(1.0))) == SVector(1.0)

    reannotated = annotate_transform(annotated; name="renamed", outputs=("v",))
    @test reannotated isa AnnotatedTransform
    @test !(reannotated.transform isa AnnotatedTransform)
    @test transform_name(reannotated) == "renamed"
    @test input_names(reannotated) == ("x",)
    @test output_names(reannotated) == ("v",)

    sequence = compose(AnnotatedTransform(Shift(1.0); inputs=("x",), outputs=("u",)), AnnotatedTransform(Scale(2.0); inputs=("u",), outputs=("v",)))
    @test input_names(sequence) == ("x",)
    @test output_names(sequence) == ("v",)
    @test sequence isa ComposedTransform

    branches = parallel(AnnotatedTransform(Shift(1.0); inputs=("x",), outputs=("u",)), AnnotatedTransform(Scale(2.0); inputs=("y",), outputs=("v",)))
    @test input_names(branches) == ("x", "y")
    @test output_names(branches) == ("u", "v")

    partially_named = parallel(AnnotatedTransform(Shift(1.0); inputs=("x",)), Scale(2.0))
    @test input_names(partially_named) === nothing
    @test output_names(partially_named) === nothing

    input = CoordinateFrame{1}(; name="input", axes_names=("x",))
    output = CoordinateFrame{1}(; name="output", axes_names=("u",))
    wcs = WCS(Step(input, annotated), Step(output))
    @test input_names(first(wcs.steps).transform) == ("x",)
    @test output_names(first(wcs.steps).transform) == ("u",)
    @test !(wcs.forward isa AnnotatedTransform)
    @test pixel_to_world(wcs, 1.0) == SVector(3.0)
    @test world_to_pixel(wcs, 3.0) == SVector(1.0)

    callable_labels = AnnotatedTransform(Shift(1.0); inputs=("model_input",), outputs=("model_output",))
    labeled_wcs = WCS(Step(input, callable_labels), Step(output))
    @test pixel_to_world(labeled_wcs, 1.0) == SVector(2.0)
    @test_throws "input names have length 2, expected 1" AnnotatedTransform(Shift(1.0); inputs=("x", "y"))
    @test_throws "output names must be unique" AnnotatedTransform(IdentityTransform{2}(); outputs=("x", "x"))
    @test_throws "transform name must not be empty" AnnotatedTransform(Shift(1.0); name="")

    compiled = compile_transform(AnnotatedTransform(AnnotationFixture.OffsetTransform(2.0); name="external", inputs=("x",), outputs=("u",)))
    @test compiled.transform == Shift(2.0)
    @test transform_name(compiled) == "external"
end
