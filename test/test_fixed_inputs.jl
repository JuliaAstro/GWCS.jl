@testset "Fixed-input transforms" begin
    affine = AffineTransform(SMatrix{2,2}(1.0, 1.0, 1.0, -1.0), SVector(0.0, 0.0))
    named = AnnotatedTransform(affine; name="sumdiff", inputs=("x", "y"), outputs=("u", "v"))
    fixed = fix_inputs(named, Dict("y" => 2.0))

    @test fixed isa FixedInputsTransform{1,2,2}
    @test transform_name(fixed) == "sumdiff"
    @test input_names(fixed) == ("x",)
    @test output_names(fixed) == ("u", "v")
    @test fixed(SVector(3.0)) == SVector(5.0, 1.0)
    @test transform_jacobian(fixed, SVector(3.0)) == SMatrix{2,1}(1.0, 1.0)

    inverse = inverse_transform(fixed)
    @test inverse(fixed(SVector(3.0))) ≈ SVector(3.0)
    @test inverse_transform(inverse) === fixed
    @test input_names(inverse) == ("u", "v")
    @test output_names(inverse) == ("x",)
    @test transform_jacobian(inverse, fixed(SVector(3.0))) ≈ SMatrix{1,2}(0.5, 0.5)
    @test all(isnan, inverse(SVector(10.0, 10.0)))

    by_index = fix_inputs(affine, 1 => 2.0)
    @test by_index(SVector(3.0)) == affine(SVector(2.0, 3.0))

    all_fixed = fix_inputs(named, (x=1.0, y=2.0))
    @test input_dimension(all_fixed) == 0
    @test input_names(all_fixed) == ()
    @test all_fixed(SVector{0,Float64}()) == named(SVector(1.0, 2.0))
    @test size(transform_jacobian(all_fixed, SVector{0,Float64}())) == (2, 0)
    all_fixed_inverse = inverse_transform(all_fixed)
    @test all_fixed_inverse(all_fixed(SVector{0,Float64}())) == SVector{0,Float64}()
    @test eltype(all_fixed_inverse(all_fixed(SVector{0,Float64}()))) === Float64

    nondual = CallableTransform{1,1}(x -> Float64(x); inverse=identity)
    fixed_nondual = fix_inputs(nondual, 1 => 2.0)
    @test size(transform_jacobian(fixed_nondual, SVector{0,Float64}())) == (1, 0)
    @test size(transform_jacobian(inverse_transform(fixed_nondual), SVector(2.0))) == (0, 1)

    no_inverse = fix_inputs(Polynomial2D([1.0 2.0; 3.0 0.0]), 2 => 1.0)
    @test !has_inverse(no_inverse)

    @test_throws "at least one input must be fixed" fix_inputs(named, Dict())
    @test_throws "same input cannot be fixed more than once" fix_inputs(named, ("x" => 1.0, 1 => 2.0))
    @test_throws "must be between 1 and 2" fix_inputs(named, 3 => 1.0)
    @test_throws "is not among" fix_inputs(named, "z" => 1.0)
    @test_throws "has no input names" fix_inputs(affine, "x" => 1.0)
    ambiguous = parallel(AnnotatedTransform(Shift(1.0); inputs=("x",)), AnnotatedTransform(Shift(2.0); inputs=("x",)))
    @test_throws "is ambiguous" fix_inputs(ambiguous, "x" => 1.0)
    @test_throws "real scalars" fix_inputs(named, "x" => [1.0])
    @test_throws "atol must be non-negative" fix_inputs(named, "x" => 1.0; atol=-1.0)

    compiled = compile_transform(fixed)
    @test compiled(SVector(3.0)) == fixed(SVector(3.0))
    @test input_names(compiled) == ("x",)

    allocated_forward(transform, coords) = @allocated transform(coords)
    allocated_inverse(transform, coords) = @allocated transform(coords)
    world = fixed(SVector(3.0))
    inverse(world)
    @test @inferred(fixed(SVector(3.0))) == world
    @test @inferred(inverse(world)) ≈ SVector(3.0)
    @test allocated_forward(fixed, SVector(3.0)) == 0
    @test allocated_inverse(inverse, world) == 0

    integer_fixed = fix_inputs(named, "y" => 2)
    integer_fixed(SVector(3.0))
    @test @inferred(integer_fixed(SVector(3.0))) == world
    @test allocated_forward(integer_fixed, SVector(3.0)) == 0
end
