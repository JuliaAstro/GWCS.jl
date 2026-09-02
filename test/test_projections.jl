function _projection_allocated(transform::T, coords::C) where {T,C}
    return @allocated transform(coords)
end

function _test_projection_performance(transform::T, coords::C) where {T,C}
    expected = transform(coords)
    @test @inferred(transform(coords)) == expected
    transform(coords)
    @test _projection_allocated(transform, coords) == 0
end

@testset "Zenithal projections" begin
    # Values generated with the corresponding Astropy projection models.
    references = (
        (GnomonicProjection(), SVector(153.43494882292202, 87.76506622828028), SVector(16.539866862653763, -28.647889756541172)),
        (StereographicProjection(), SVector(153.43494882292202, 87.7642157682017), SVector(15.352357850242075, -26.591063812598186)),
        (ZenithalEquidistant(), SVector(153.43494882292202, 87.76393202250021), SVector(15.0, -25.98076211353316)),
        (ZenithalEqualArea(), SVector(153.43494882292202, 87.76379009291757), SVector(14.829238941980536, -25.684995285089236)),
        (SlantOrthographic(0.1, 0.2), SVector(153.43494882292202, 87.7730474413347), SVector(15.091562770782687, -23.27456450895644)),
        (ZenithalPerspective(1.5, 15.0), SVector(152.63219484137733, 87.8326085109999), SVector(14.427529882336112, -25.870738832948405)),
        (SlantZenithalPerspective(1.5, 15.0, 80.0), SVector(153.5130711980779, 87.7606067237832), SVector(15.36376104106826, -27.052503716952117)),
        (AiryProjection(75.8), SVector(153.43494882292202, 87.7553857871437), SVector(15.030376707991271, -26.03337611514073))
    )

    intermediate = SVector(1.0, 2.0)
    native = SVector(30.0, 60.0)
    for (projection, expected_native, expected_intermediate) in references
        to_native = IntermediateToNative(projection)
        to_intermediate = NativeToIntermediate(projection)
        @test to_native(intermediate) ≈ expected_native atol=5e-13
        @test to_intermediate(native) ≈ expected_intermediate atol=5e-13
        @test to_intermediate(to_native(intermediate)) ≈ intermediate atol=1e-10
        @test to_native(to_intermediate(native)) ≈ native atol=1e-10
        @test inverse_transform(to_native) == to_intermediate
        @test inverse_transform(to_intermediate) == to_native
        @test all(axis_correlation_matrix(to_native))
        _test_projection_performance(to_native, intermediate)
        _test_projection_performance(to_intermediate, native)
    end

    @test IntermediateToNative{GnomonicProjection}()(intermediate) == IntermediateToNative(GnomonicProjection())(intermediate)
    @test ZenithalPerspective{Float64}(1, 0).mu === 1.0
    @test SlantZenithalPerspective{Float64}(1, 0, 90).theta0 === 90.0
    @test SlantOrthographic{Float64}(1, 0).xi === 1.0
    @test AiryProjection{Float64}(90).theta_b === 90.0

    float32_projections = (
        GnomonicProjection(),
        StereographicProjection(),
        ZenithalEquidistant(),
        ZenithalEqualArea(),
        SlantOrthographic(0.1f0, 0.2f0),
        ZenithalPerspective(1.5f0, 15.0f0),
        SlantZenithalPerspective(1.5f0, 15.0f0, 80.0f0),
        AiryProjection(75.8f0)
    )
    for projection in float32_projections
        @test eltype(IntermediateToNative(projection)(SVector(1.0f0, 2.0f0))) === Float32
        @test eltype(NativeToIntermediate(projection)(SVector(30.0f0, 60.0f0))) === Float32
    end

    @test_throws "undefined for mu = -1" ZenithalPerspective(-1)
    @test_throws "undefined when cos(gamma) = 0" ZenithalPerspective(0, 90)
    @test_throws "requires mu * sin(theta0) + 1 != 0" SlantZenithalPerspective(-1, 0, 90)
    @test_throws "greater than -90" AiryProjection(-90)
    @test_throws DomainError NativeToIntermediate(GnomonicProjection())(SVector(0.0, 0.0))
    @test_throws DomainError NativeToIntermediate(StereographicProjection())(SVector(0.0, -90.0))
    @test_throws DomainError IntermediateToNative(ZenithalEquidistant())(SVector(0.0, 181.0))
    @test_throws DomainError IntermediateToNative(ZenithalEqualArea())(SVector(0.0, 120.0))
    @test_throws DomainError IntermediateToNative(SlantOrthographic())(SVector(60.0, 0.0))
    @test_throws DomainError NativeToIntermediate(SlantOrthographic())(SVector(0.0, -1.0))
    @test_throws DomainError IntermediateToNative(ZenithalPerspective(1.5, 15.0))(SVector(-500.0, -500.0))
    @test_throws DomainError NativeToIntermediate(ZenithalPerspective(1.5, 15.0))(SVector(-180.0, -90.0))
    @test_throws DomainError IntermediateToNative(SlantZenithalPerspective(1.5, 15.0, 80.0))(SVector(-500.0, -500.0))
    @test_throws DomainError NativeToIntermediate(SlantZenithalPerspective(1.5, 15.0, 80.0))(SVector(-180.0, -90.0))
    @test_throws DomainError NativeToIntermediate(AiryProjection())(SVector(0.0, -90.0))

    airy = IntermediateToNative(AiryProjection(75.8))
    finite_difference = (airy(SVector(1.0, 2.0 + 1e-6)) - airy(SVector(1.0, 2.0 - 1e-6))) / 2e-6
    @test transform_jacobian(airy, SVector(1.0, 2.0))[:, 2] ≈ finite_difference rtol=1e-6
end

@testset "Celestial projection construction" begin
    @test default_lonpole(20.0, 0.0) == 0.0
    @test default_lonpole(-20.0, 0.0) == 180.0
    @test default_lonpole(-5.0, 90.0) == 180.0
    @test SVector(GWCS._native_pole(30.0, 20.0, 0.0, 0.0, 0.0, 90.0)) ≈ SVector(210.0, 70.0)
    @test SVector(GWCS._native_pole(30.0, -20.0, 0.0, 0.0, 180.0, 90.0)) ≈ SVector(30.0, 70.0)
    @test_throws "no valid native-pole latitude exists" GWCS._native_pole(30.0, 20.0, 0.0, 0.0, 120.0, -30.0)

    transform = celestial_transform(GnomonicProjection(), (83.8221, -5.3911))
    @test transform(SVector(0.0, 0.0)) ≈ SVector(83.8221, -5.3911) atol=1e-12
    @test inverse_transform(transform)(transform(SVector(1.0, 2.0))) ≈ SVector(1.0, 2.0) atol=1e-10
    @test celestial_transform(ZenithalEqualArea(), [10.0, 20.0]; lonpole=0.0)(SVector(0.0, 0.0)) ≈ SVector(10.0, 20.0) atol=1e-12
    @test_throws DimensionMismatch celestial_transform(GnomonicProjection(), (1.0,))
    @test_throws "must be real numbers" celestial_transform(GnomonicProjection(), (0.0, "north"))
    @test_throws "lonpole must be a real number or nothing" celestial_transform(GnomonicProjection(), (0.0, 0.0); lonpole="north")
    @test_throws "reference latitude must be between -90 and 90 degrees" celestial_transform(GnomonicProjection(), (0.0, 91.0))
    _test_projection_performance(transform, SVector(1.0, 2.0))
    _test_projection_performance(inverse_transform(transform), transform(SVector(1.0, 2.0)))

    azp_transform = celestial_transform(ZenithalPerspective(1.5, 15.0), (83.8221, -5.3911))
    wcs = WCS(Step(Frame2D(; name="intermediate"), azp_transform), Step(CelestialFrame(; name="world")))
    world = pixel_to_world(wcs, SVector(1.0, 2.0))
    @test world_to_pixel(wcs, world) ≈ SVector(1.0, 2.0) atol=1e-10
    @test @inferred(pixel_to_world(wcs, SVector(1.0, 2.0))) == world
end
