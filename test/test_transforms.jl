@testset "Leaf transforms" begin
    @test input_dimension(Shift(2.0)) == 1
    @test output_dimension(Shift(2.0)) == 1
    @test Shift(2.0)(SVector(3.0)) == SVector(5.0)
    @test inverse_transform(Shift(2.0))(SVector(5.0)) == SVector(3.0)
    @test Scale(3.0)(SVector(2.0)) == SVector(6.0)
    @test inverse_transform(Scale(3.0))(SVector(6.0)) == SVector(2.0)
    @test !has_inverse(Scale(0.0))
    @test eltype(Shift(2.0f0)(SVector(3.0f0))) === Float32

    constant = ConstantTransform{2}(42.0)
    @test input_dimension(constant) == 2
    @test output_dimension(constant) == 1
    @test constant(SVector(1.0, 2.0)) == SVector(42.0)
    @test transform_jacobian(constant, SVector(1.0, 2.0)) == SMatrix{1,2}(0.0, 0.0)
    @test !has_inverse(constant)

    inserted = ConstantTransform(3.0)
    @test input_dimension(inserted) == 0
    @test inserted(SVector{0,Float64}()) == SVector(3.0)
    @test size(transform_jacobian(inserted, SVector{0,Float64}())) == (1, 0)
    inserted_inverse = inverse_transform(inserted)
    @test inserted_inverse(SVector(3.0)) == SVector{0,Float64}()
    @test inverse_transform(inserted_inverse) === inserted
    @test_throws DomainError inserted_inverse(SVector(4.0))

    add_axis = parallel(IdentityTransform{2}(), inserted)
    @test input_dimension(add_axis) == 2
    @test output_dimension(add_axis) == 3
    @test add_axis(SVector(1.0, 2.0)) == SVector(1.0, 2.0, 3.0)
    @test inverse_transform(add_axis)(SVector(1.0, 2.0, 3.0)) == SVector(1.0, 2.0)
    @test_throws DomainError inverse_transform(add_axis)(SVector(1.0, 2.0, 4.0))

    affine = AffineTransform(SMatrix{2,2}(2.0, -1.0, 1.0, 3.0), SVector(4.0, -2.0))
    coords = SVector(1.0, 2.0)
    result = affine(coords)
    @test result == SVector(8.0, 3.0)
    @test inverse_transform(affine)(result) ≈ coords

    remap = RemapAxes{2}((1, 2, 1, 2))
    @test remap(coords) == SVector(1.0, 2.0, 1.0, 2.0)
    @test inverse_transform(remap)(remap(coords)) == coords
    @test_throws "mapping indices must be between 1 and 2" RemapAxes{2}((1, 3))

    poly1 = Polynomial1D([1.0, 2.0, 3.0])
    @test degree(poly1) == 2
    @test poly1(SVector(2.0)) == SVector(17.0)

    coefficients = [1.0 3.0 6.0; 2.0 5.0 0.0; 4.0 0.0 0.0]
    poly2 = Polynomial2D(coefficients)
    @test degree(poly2) == 2
    @test poly2(SVector(2.0, 3.0)) == SVector(1 + 2 * 2 + 4 * 2^2 + 3 * 3 + 5 * 2 * 3 + 6 * 3^2)
    @test ForwardDiff.jacobian(x -> poly2(SVector{2}(x)), [2.0, 3.0]) ≈ [33.0 49.0]

    spherical = SphericalToCartesian(360)
    cartesian = spherical(SVector(90.0, 0.0))
    @test cartesian ≈ SVector(0.0, 1.0, 0.0) atol=1e-15
    @test inverse_transform(spherical)(cartesian) ≈ SVector(90.0, 0.0) atol=1e-14

    tangent = TangentPlaneToNative()
    native = tangent(SVector(1.0, 2.0))
    @test inverse_transform(tangent)(native) ≈ SVector(1.0, 2.0) atol=1e-14

    celestial = NativeToCelestial(83.8221, -5.3911, 180.0)
    native_coords = SVector(0.0, 90.0)
    @test celestial(native_coords) ≈ SVector(83.8221, -5.3911) atol=1e-12
    @test inverse_transform(celestial)(celestial(SVector(10.0, 80.0))) ≈ SVector(10.0, 80.0) atol=1e-12

    rotation = Rotation3D([90.0], "z")
    @test rotation(SVector(1.0, 0.0, 0.0)) ≈ SVector(0.0, -1.0, 0.0) atol=1e-15
    coords3 = SVector(0.2, 0.3, 0.4)
    @test inverse_transform(rotation)(rotation(coords3)) ≈ coords3

    explicit = WithInverse(poly1, Polynomial1D([-1 / 3, 1 / 3]))
    @test has_inverse(explicit)
    @test explicit(SVector(2.0)) == SVector(17.0)
end
