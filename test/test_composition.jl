struct ArithmeticUncompiled <: AbstractTransform{1,1} end
struct ArithmeticCompiled <: AbstractTransform{1,1} end

@inline (::ArithmeticUncompiled)(coords::StaticVector{1}) = coords
@inline (::ArithmeticCompiled)(coords::StaticVector{1}) = coords

GWCS.compile_transform(::ArithmeticUncompiled) = ArithmeticCompiled()

@testset "Transform composition" begin
    sequence = compose(Shift(2.0), Scale(3.0), Shift(-1.0))
    @test sequence(SVector(4.0)) == SVector(17.0)
    @test inverse_transform(sequence)(sequence(SVector(4.0))) ≈ SVector(4.0)
    @test (Shift(2.0) | Scale(3.0))(SVector(4.0)) == SVector(18.0)

    affine_parallel = parallel(Shift(1.0), Scale(2.0))
    @test affine_parallel isa AffineTransform
    @test affine_parallel(SVector(3.0, 4.0)) == SVector(4.0, 8.0)
    @test (Shift(1.0) & Scale(2.0))(SVector(3.0, 4.0)) == SVector(4.0, 8.0)

    polynomial = Polynomial2D([1.0 2.0; 3.0 0.0])
    polys = parallel(polynomial, polynomial)
    @test polys(SVector(2.0, 3.0, 4.0, 5.0)) == SVector(13.0, 23.0)

    remap = RemapAxes{2}((1, 2, 1, 2))
    pair = parallel(polynomial, polynomial)
    mapped = compose(remap, pair)
    @test mapped(SVector(2.0, 3.0)) == SVector(13.0, 13.0)
    @test !has_inverse(mapped)

    @test_throws DimensionMismatch compose(Shift(1.0), Rotation3D([1.0], "x"))
    @test_throws DimensionMismatch compose(Shift(1.0), IdentityTransform{2}())
    @test_throws DimensionMismatch compose(IdentityTransform{2}(), Shift(1.0))
    @test compose(IdentityTransform{2}()) isa IdentityTransform{2}
end

@testset "Elementwise arithmetic composition" begin
    left = Shift(1.0)
    right = Scale(2.0)
    coordinate = SVector(3.0)

    @test elementwise(+, left, right) isa ArithmeticTransform
    @test (left + right)(coordinate) == SVector(10.0)
    @test (left - right)(coordinate) == SVector(-2.0)
    @test (left * right)(coordinate) == SVector(24.0)
    @test (left / right)(coordinate) == SVector(2 / 3)
    @test (left ^ right)(coordinate) == SVector(4096.0)
    @test isinf(only((IdentityTransform{1}() / ConstantTransform{1}(0.0))(SVector(1.0))))
    @test_throws DomainError (Shift(-2.0) ^ ConstantTransform{1}(0.5))(SVector(1.0))

    multi_left = AffineTransform((@SMatrix [1.0 2.0; 3.0 4.0]), (@SVector [1.0; -1.0]))
    multi_right = AffineTransform(@SMatrix [2.0 0.0; 0.0 0.5])
    @test (multi_left + multi_right)(SVector(2.0, 3.0)) == SVector(13.0, 18.5)
    @test (IdentityTransform{2}() * multi_right)(SVector(2.0, 3.0)) == SVector(8.0, 4.5)

    integer_sum = IdentityTransform{1}() + IdentityTransform{1}()
    @test integer_sum(SVector(2)) == SVector(4)
    @test (left + right)(SVector(3 // 1)) == SVector(10.0)

    @test_throws "unsupported elementwise transform operator" elementwise(max, left, right)
    @test_throws "input dimensions 1 and 2" left + IdentityTransform{2}()
    @test_throws "output dimensions 1 and 2" Polynomial2D([1.0 0.0; 0.0 0.0]) + IdentityTransform{2}()
    @test !has_inverse(left + right)

    named_left = AnnotatedTransform(left; inputs=("x",), outputs=("y",))
    named_right = AnnotatedTransform(right; inputs=("x",), outputs=("y",))
    unnamed = left + right
    @test input_names(named_left + named_right) == ("x",)
    @test output_names(named_left + unnamed) == ("y",)
    @test input_names(named_left + AnnotatedTransform(right; inputs=("u",))) === nothing
    @test output_names(named_left + AnnotatedTransform(right; outputs=("v",))) === nothing

    compiled = compile_transform(ArithmeticUncompiled() + ArithmeticUncompiled())
    @test compiled.left isa ArithmeticCompiled
    @test compiled.right isa ArithmeticCompiled

    for transform in (left + right, left - right, left * right, left / right, left ^ right)
        @test transform_jacobian(transform, coordinate) ≈ ForwardDiff.jacobian(transform, coordinate)
    end
    negative_square = Shift(-3.0) ^ ConstantTransform{1}(2.0)
    @test negative_square(SVector(1.0)) == SVector(4.0)
    @test transform_jacobian(negative_square, SVector(1.0)) == SMatrix{1,1}(-4.0)
    zero_square = Shift(-1.0) ^ ConstantTransform{1}(2.0)
    @test transform_jacobian(zero_square, SVector(1.0)) == SMatrix{1,1}(0.0)

    arithmetic = left * right
    allocated_arithmetic(transform, value) = @allocated transform(value)
    arithmetic(coordinate)
    @test @inferred(arithmetic(coordinate)) == SVector(24.0)
    @test allocated_arithmetic(arithmetic, coordinate) == 0
end
