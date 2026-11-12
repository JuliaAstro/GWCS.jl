using GWCS: AbstractTransform, AffineTransform, AnnotatedTransform, CallableTransform, CartesianToSpherical, ConstantTransform, FITSLookupDistortion,
            FITSTabularTable, FITSTabularTransform, GnomonicProjection, IdentityTransform, IntermediateToNative, Polynomial1D, Polynomial2D,
            RemapAxes, Rotation3D, Scale, Shift, SphericalToCartesian, WithInverse, compose, fix_inputs, inverse_transform, parallel

struct DependencyFallback <: AbstractTransform{2,1} end
struct DependencyOverride <: AbstractTransform{2,2} end
struct BadDependencyShape <: AbstractTransform{2,1} end
struct BadDependencyValues <: AbstractTransform{1,1} end

@inline (::DependencyFallback)(coords::StaticVector{2}) = SVector(coords[1] + coords[2])
@inline (::DependencyOverride)(coords::StaticVector{2}) = coords
@inline (::BadDependencyShape)(coords::StaticVector{2}) = SVector(coords[1])
@inline (::BadDependencyValues)(coords::StaticVector{1}) = coords

GWCS.transform_dependencies(::DependencyOverride) = @SMatrix [true false; false true]
GWCS.transform_dependencies(::BadDependencyShape) = @SMatrix [true]
GWCS.transform_dependencies(::BadDependencyValues) = SMatrix{1,1}(2)

@testset "Transform dependencies" begin
    @test axis_correlation_matrix(IdentityTransform{2}()) == @SMatrix [true false; false true]
    @test axis_correlation_matrix(Shift(1.0)) == @SMatrix [true]
    @test axis_correlation_matrix(Scale(0.0)) == @SMatrix [false]
    @test axis_correlation_matrix(Scale(2.0)) == @SMatrix [true]
    @test axis_correlation_matrix(AffineTransform(@SMatrix [1.0 0.0; 2.0 0.0])) == @SMatrix [true false; true false]
    @test axis_correlation_matrix(RemapAxes{3}((3, 1))) == @SMatrix [false false true; true false false]
    @test axis_correlation_matrix(ConstantTransform{2}(1.0)) == @SMatrix [false false]
    @test size(axis_correlation_matrix(inverse_transform(ConstantTransform(1.0)))) == (0, 1)

    @test axis_correlation_matrix(Polynomial1D([1.0])) == @SMatrix [false]
    @test axis_correlation_matrix(Polynomial1D([1.0, 2.0])) == @SMatrix [true]
    @test axis_correlation_matrix(Polynomial2D([1.0 2.0; 0.0 0.0])) == @SMatrix [false true]
    @test axis_correlation_matrix(Polynomial2D([1.0 0.0; 2.0 0.0])) == @SMatrix [true false]
    @test axis_correlation_matrix(Polynomial2D([1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 0.0])) == @SMatrix [true true]

    @test axis_correlation_matrix(SphericalToCartesian()) == @SMatrix [true true; true true; false true]
    @test axis_correlation_matrix(CartesianToSpherical()) == @SMatrix [true true false; true true true]
    @test all(axis_correlation_matrix(IntermediateToNative(GnomonicProjection())))
    @test axis_correlation_matrix(Rotation3D(@SMatrix [1.0 0.0 0.0; 0.0 0.0 -1.0; 0.0 1.0 0.0])) == @SMatrix [true false false; false false true; false true false]

    composed = compose(RemapAxes{2}((2, 1)), AffineTransform(@SMatrix [1.0 0.0; 0.0 0.0]))
    @test axis_correlation_matrix(composed) == @SMatrix [false true; false false]
    block = parallel(Polynomial1D([0.0, 1.0]), Polynomial1D([1.0]))
    @test axis_correlation_matrix(block) == @SMatrix [true false; false false]
    arithmetic_left = AffineTransform(@SMatrix [1.0 0.0; 0.0 0.0])
    arithmetic_right = AffineTransform(@SMatrix [0.0 0.0; 0.0 1.0])
    @test axis_correlation_matrix(arithmetic_left + arithmetic_right) == @SMatrix [true false; false true]

    affine = AffineTransform(@SMatrix [1.0 1.0; 1.0 -1.0])
    fixed = fix_inputs(affine, 2 => 1.0)
    @test axis_correlation_matrix(fixed) == @SMatrix [true; true]
    @test axis_correlation_matrix(inverse_transform(fixed)) == @SMatrix [true true]
    @test axis_correlation_matrix(AnnotatedTransform(affine; inputs=("x", "y"))) == @SMatrix [true true; true true]
    paired = WithInverse(Scale(0.0), Shift(1.0))
    @test axis_correlation_matrix(paired) == @SMatrix [false]
    @test axis_correlation_matrix(inverse_transform(paired)) == @SMatrix [true]

    declared = CallableTransform{2,2}((x, y) -> (x, x + y); inverse=(u, v) -> (u, v - u),
        dependencies=(@SMatrix [true false; true true]), inverse_dependencies=(@SMatrix [true false; true true]))
    @test axis_correlation_matrix(declared) == @SMatrix [true false; true true]
    @test axis_correlation_matrix(inverse_transform(declared)) == @SMatrix [true false; true true]
    @test all(axis_correlation_matrix(CallableTransform{2,2}((x, y) -> (x, y))))
    @test_throws "has size (1, 2), expected (2, 2)" CallableTransform{2,2}((x, y) -> (x, y); dependencies=(@SMatrix [true false]))
    @test_throws "requires an inverse callable" CallableTransform{1,1}(identity; inverse_dependencies=(@SMatrix [true]))

    @test all(axis_correlation_matrix(DependencyFallback()))
    @test axis_correlation_matrix(DependencyOverride()) == @SMatrix [true false; false true]
    @test_throws DimensionMismatch axis_correlation_matrix(BadDependencyShape())
    @test_throws "entries must be boolean-convertible" axis_correlation_matrix(BadDependencyValues())

    # FITS lookup and TAB transforms depend only on the axes their tables read.
    table = GWCS.FITSLookupTable([0.0, 1.0, 2.0], (1,), Val(3))
    @test axis_correlation_matrix(table) == SMatrix{1,3,Bool}(true, false, false)
    lookup = FITSLookupDistortion(Val(3), (table, nothing), IdentityTransform{3}(), (nothing, table))
    @test axis_correlation_matrix(lookup) == SMatrix{3,3,Bool}([1 0 0; 1 1 0; 0 0 1])
    tabular = FITSTabularTransform(Val(3), (FITSTabularTable((1, 3), reshape(collect(1.0:8.0), 2, 2, 2), ([1.0, 2.0], [1.0, 2.0])),), (1.0, 1.0, 1.0))
    @test axis_correlation_matrix(tabular) == SMatrix{3,3,Bool}([1 0 1; 0 1 0; 1 0 1])
    separable = FITSTabularTransform(Val(2), (FITSTabularTable((2,), reshape([1.0, 2.0, 4.0], 1, 3), ([1.0, 2.0, 3.0],)),), (1.0, 1.0))
    @test axis_correlation_matrix(inverse_transform(separable)) == SMatrix{2,2,Bool}([1 0; 0 1])
end
