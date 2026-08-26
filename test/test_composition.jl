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
