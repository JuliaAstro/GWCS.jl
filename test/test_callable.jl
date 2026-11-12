using GWCS: CallableTransform, CoordinateFrame, Shift, Step, has_inverse, input_dimension, inverse_transform, output_dimension, transform,
            transform_jacobian

module CallableFixture

using ForwardDiff

forward(x, y) = (x + y, x - y)
backward(u, v) = ((u + v) / 2, (u - v) / 2)
jacobian(x, y) = [1 1; 1 -1]

function nonlinear(x)
    x isa ForwardDiff.Dual && error("analytical Jacobian was not used")
    return x + 0.1 * x^2
end

nonlinear_jacobian(x) = 1 + 0.2 * x

end

@testset "Callable transforms" begin
    transform = CallableTransform{2,2}(CallableFixture.forward; inverse=CallableFixture.backward, jacobian=CallableFixture.jacobian)
    coords = SVector(3.0, 1.0)
    outputs = transform(coords)
    @test input_dimension(transform) == 2
    @test output_dimension(transform) == 2
    @test outputs == SVector(4.0, 2.0)
    @test inverse_transform(transform)(outputs) == coords
    @test transform_jacobian(transform, coords) == SMatrix{2,2}(1, 1, 1, -1)
    @test transform_jacobian(inverse_transform(transform), outputs) == SMatrix{2,2}(0.5, 0.5, 0.5, -0.5)
    @test inverse_transform(inverse_transform(transform)) === transform
    @test transform_jacobian(inverse_transform(inverse_transform(transform)), coords) == SMatrix{2,2}(1, 1, 1, -1)

    outer = CallableTransform(CallableFixture.forward, 2, 2)
    @test outer(coords) == outputs
    @test_throws "input dimension must be non-negative" CallableTransform(identity, -1, 1)
    @test_throws "output dimension must be positive" CallableTransform(identity, 1, 0)

    offset = 2.5
    closure = CallableTransform{1,1}(x -> x + offset)
    @test closure(SVector(1.0)) == SVector(3.5)
    @test !has_inverse(closure)
    @test transform_jacobian(closure, SVector(1.0)) == SMatrix{1,1}(1.0)

    float32_transform = CallableTransform{1,1}(x -> x * 2)
    @test eltype(float32_transform(SVector(1.5f0))) === Float32
    @test ForwardDiff.derivative(x -> only(float32_transform(SVector(x))), 1.5) == 2

    vector_output = CallableTransform{2,2}((x, y) -> [x, y])
    @test vector_output(coords) == coords
    @test_throws "returned one scalar, expected 2 outputs" CallableTransform{1,2}(identity)(SVector(1.0))
    @test_throws "returned 1 outputs, expected 2" CallableTransform{1,2}(x -> (x,))(SVector(1.0))
    @test_throws "expected a scalar, tuple, or AbstractVector" CallableTransform{1,1}(x -> "bad")(SVector(1.0))

    scalar_jacobian = CallableTransform{1,1}(x -> x^2; jacobian=x -> 2x)
    @test transform_jacobian(scalar_jacobian, SVector(3.0)) == SMatrix{1,1}(6.0)
    @test_throws "expected an AbstractMatrix" transform_jacobian(CallableTransform{1,2}(x -> (x, x); jacobian=x -> (1, 1)), SVector(1.0))
    @test_throws "has size (1, 2), expected (2, 1)" transform_jacobian(CallableTransform{1,2}(x -> (x, x); jacobian=x -> [1 1]), SVector(1.0))

    first = CallableTransform{2,2}(CallableFixture.forward; jacobian=CallableFixture.jacobian)
    second = CallableTransform{2,2}((u, v) -> (2u, 3v); jacobian=(u, v) -> [2 0; 0 3])
    @test transform_jacobian(first | second, coords) == SMatrix{2,2}(2, 3, 2, -3)

    parallel_transform = CallableTransform{1,1}(x -> x^2; jacobian=x -> 2x) & CallableTransform{1,1}(x -> 3x; jacobian=x -> 3)
    @test transform_jacobian(parallel_transform, SVector(2.0, 4.0)) == SMatrix{2,2}(4.0, 0.0, 0.0, 3.0)
    allocated_jacobian(transform, value) = @allocated transform_jacobian(transform, value)
    transform_jacobian(parallel_transform, SVector(2.0, 4.0))
    @test allocated_jacobian(parallel_transform, SVector(2.0, 4.0)) == 0

    pixels = CoordinateFrame{1}(name="pixels")
    world = CoordinateFrame{1}(name="world")
    nonlinear = CallableTransform{1,1}(CallableFixture.nonlinear; jacobian=CallableFixture.nonlinear_jacobian)
    nonlinear_wcs = WCS(Step(pixels, nonlinear | Shift(2.0)), Step(world))
    pixel = SVector(1.5)
    transformed = pixel_to_world(nonlinear_wcs, pixel)
    @test world_to_pixel(nonlinear_wcs, transformed; tolerance=1e-12, quiet=false) ≈ pixel atol=1e-12

    compiled_wcs = WCS(Step(pixels, CallableTransform{1,1}(x -> x + 2; inverse=x -> x - 2)), Step(world))
    allocated_forward(wcs, value) = @allocated pixel_to_world(wcs, value)
    allocated_inverse(wcs, value) = @allocated world_to_pixel(wcs, value)
    pixel_to_world(compiled_wcs, SVector(1.0))
    world_to_pixel(compiled_wcs, SVector(3.0))
    @test @inferred(pixel_to_world(compiled_wcs, SVector(1.0))) == SVector(3.0)
    @test @inferred(world_to_pixel(compiled_wcs, SVector(3.0))) == SVector(1.0)
    @test allocated_forward(compiled_wcs, SVector(1.0)) == 0
    @test allocated_inverse(compiled_wcs, SVector(3.0)) == 0
end
