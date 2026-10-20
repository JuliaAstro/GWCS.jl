function _projection_allocated(transform::T, coords::C) where {T,C}
    return @allocated transform(coords)
end

function _projection_jacobian_allocated(transform::T, coords::C) where {T,C}
    return @allocated transform_jacobian(transform, coords)
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
    @test all(isnan, NativeToIntermediate(GnomonicProjection())(SVector(0.0, 0.0)))
    @test all(isnan, NativeToIntermediate(StereographicProjection())(SVector(0.0, -90.0)))
    @test all(isnan, IntermediateToNative(ZenithalEquidistant())(SVector(0.0, 181.0)))
    @test all(isnan, IntermediateToNative(ZenithalEqualArea())(SVector(0.0, 120.0)))
    @test all(isnan, IntermediateToNative(SlantOrthographic())(SVector(60.0, 0.0)))
    @test all(isnan, NativeToIntermediate(SlantOrthographic())(SVector(0.0, -1.0)))
    @test all(isnan, IntermediateToNative(ZenithalPerspective(1.5, 15.0))(SVector(-500.0, -500.0)))
    @test all(isnan, NativeToIntermediate(ZenithalPerspective(1.5, 15.0))(SVector(-180.0, -90.0)))
    @test all(isnan, IntermediateToNative(SlantZenithalPerspective(1.5, 15.0, 80.0))(SVector(-500.0, -500.0)))
    @test all(isnan, NativeToIntermediate(SlantZenithalPerspective(1.5, 15.0, 80.0))(SVector(-180.0, -90.0)))
    @test all(isnan, NativeToIntermediate(AiryProjection())(SVector(0.0, -90.0)))

    cyp = CylindricalPerspective(2.0, 3.0)
    @test cyp.mu == 2.0
    @test cyp.lambda == 3.0
    @test IntermediateToNative(cyp)(SVector(1.0, 2.0)) ≈ SVector(0.33333333333333331, 1.2000000012828205) atol=2e-15
    @test NativeToIntermediate(cyp)(SVector(30.0, 60.0)) ≈ SVector(90.0, 99.23920117592256) atol=2e-14
    @test_throws "lambda must be nonzero" CylindricalPerspective(1, 0)

    airy = IntermediateToNative(AiryProjection(75.8))
    finite_difference = (airy(SVector(1.0, 2.0 + 1e-6)) - airy(SVector(1.0, 2.0 - 1e-6))) / 2e-6
    @test transform_jacobian(airy, SVector(1.0, 2.0))[:, 2] ≈ finite_difference rtol=1e-6

    zpn = IntermediateToNative(ZenithalPolynomial((0.0, 1.0, 0.05, 0.01)))
    finite_difference = (zpn(SVector(1.0, 2.0 + 1e-6)) - zpn(SVector(1.0, 2.0 - 1e-6))) / 2e-6
    @test transform_jacobian(zpn, SVector(1.0, 2.0))[:, 2] ≈ finite_difference rtol=1e-6

    pco_to_native = IntermediateToNative(PolyconicProjection())
    pco_to_intermediate = NativeToIntermediate(PolyconicProjection())
    @test transform_jacobian(pco_to_native, SVector(30.0, 0.0)) ≈ (@SMatrix [1.0 0.0; 0.0 inv(1 + deg2rad(30.0)^2 / 2)]) atol=2e-15
    @test transform_jacobian(pco_to_native, SVector(0.0, -89.0)) ≈ (@SMatrix [inv(cosd(89.0)) 0.0; 0.0 1.0]) atol=2e-13
    @test transform_jacobian(pco_to_intermediate, SVector(-90.0, -0.01)) ≈ (@SMatrix [0.9999999471884624 -0.0004996396702968213; 0.00027415566880621375 2.233700451793191]) atol=2e-15
    pco_jacobian = @inferred transform_jacobian(pco_to_native, SVector(1.0f0, 2.0f0))
    @test eltype(pco_jacobian) === Float32
    transform_jacobian(pco_to_native, SVector(1.0, 2.0))
    @test _projection_jacobian_allocated(pco_to_native, SVector(1.0, 2.0)) == 0
end

@testset "WCSLIB projection boundaries" begin
    # Sparse boundary/extreme cases generated directly with WCSLIB 8.9 prjs2x/prjx2s.
    cube_references = (
        (TangentialSphericalCube(), SVector(45.0, 0.0), SVector(44.999999999999993, 0.0)),
        (TangentialSphericalCube(), SVector(0.0, 45.0), SVector(0.0, 44.999999999999993)),
        (COBESphericalCube(), SVector(45.0, 0.0), SVector(45.0, 0.0)),
        (QuadrilateralizedSphericalCube(), SVector(45.0, 0.0), SVector(44.999999999999993, 0.0)),
    )
    for (projection, native, expected) in cube_references
        forward = NativeToIntermediate(projection)
        @test forward(native) ≈ expected atol=2e-14
        @test IntermediateToNative(projection)(forward(native)) ≈ native atol=2e-14
    end

    quadratic_zpn = ZenithalPolynomial((0.0, 1.0, -0.2))
    @test quadratic_zpn.monotonic_limit == 2.5
    @test quadratic_zpn.radial_limit == 1.25
    quadratic_boundary = SVector(0.0, 90.0 - rad2deg(2.5 * (1 - 1e-6)))
    @test IntermediateToNative(quadratic_zpn)(NativeToIntermediate(quadratic_zpn)(quadratic_boundary)) ≈ quadratic_boundary atol=3e-7
    @test all(isnan, NativeToIntermediate(quadratic_zpn)(SVector(0.0, 90.0 - rad2deg(2.5) - 0.01)))
    @test all(isnan, IntermediateToNative(quadratic_zpn)(SVector(0.0, rad2deg(1.25) + 0.01)))

    strong_zpn = ZenithalPolynomial((0.0, 0.95, -0.025, -0.15833, 0.00208, 0.00792, -0.00007, -0.00019))
    extreme_references = (
        (strong_zpn, SVector(179.0, 10.0), SVector(0.89412934323739812, 51.22463576700688)),
        (ZenithalPerspective(0.5, 30.0), SVector(170.0, 5.0), SVector(716.82402686439639, 4694.2168828492022)),
        (ConicPerspective(60.0, 15.0), SVector(170.0, -24.999999), SVector(359.74490843286873, 590.68807277234828)),
    )
    for (projection, native, expected) in extreme_references
        projected = NativeToIntermediate(projection)(native)
        @test projected ≈ expected rtol=3e-14
        @test IntermediateToNative(projection)(projected) ≈ native atol=3e-11
    end
    zpn_limit = 90.0 - rad2deg(strong_zpn.monotonic_limit)
    zpn_boundary = SVector(0.0, 90.0 - rad2deg(strong_zpn.monotonic_limit * (1 - 1e-6)))
    projected_boundary = NativeToIntermediate(strong_zpn)(zpn_boundary)
    @test IntermediateToNative(strong_zpn)(projected_boundary) ≈ zpn_boundary atol=2e-7
    @test all(isnan, NativeToIntermediate(strong_zpn)(SVector(0.0, zpn_limit - 0.01)))
    @test all(isnan, IntermediateToNative(strong_zpn)(SVector(0.0, rad2deg(strong_zpn.radial_limit) * (1 + 1e-6))))

    car = NativeToIntermediate(PlateCarreeProjection())
    @test car(SVector(180.0, 23.0)) == SVector(180.0, 23.0)
    @test car(SVector(-180.0, 23.0)) == SVector(-180.0, 23.0)
    @test car(SVector(181.0, 23.0)) ≈ SVector(-179.0, 23.0) atol=3e-14
    @test car(SVector(540.0, 23.0)) == SVector(-180.0, 23.0)
end

@testset "Projection parameters and non-finite coordinates" begin
    # Conic and Bonne constants are precomputed at construction, so invalid parameters fail there.
    @test ConicEquidistant{Float32}(45, 10) isa ConicEquidistant{Float32}
    @test BonneProjection{Float32}(30).theta1 === 30.0f0
    @test_throws "COP is degenerate for sigma = 0" ConicPerspective(0.0, 0.0)
    @test_throws "COP is degenerate for delta" ConicPerspective(45.0, 90.0)
    @test_throws "COD conic constant is degenerate" ConicEquidistant(0.0, 0.0)
    @test_throws "between -90 and 90 degrees" ConicOrthomorphic(91.0, 0.0)
    @test_throws "Bonne standard parallel" BonneProjection(91.0)
    for projection in (ConicPerspective(60.0, 15.0), ConicEquidistant(45.0, 25.0), ConicEqualArea(-45.0, 10.0), ConicOrthomorphic(45.0, 10.0), BonneProjection(-45.0))
        native = SVector(20.0, 30.0)
        _test_projection_performance(NativeToIntermediate(projection), native)
        _test_projection_performance(IntermediateToNative(projection), NativeToIntermediate(projection)(native))
    end

    for projection in (GnomonicProjection(), PlateCarreeProjection(), HammerAitoffProjection(), MollweideProjection(), ZenithalPolynomial((0.0, 1.0, -0.2)), ConicEqualArea(), HEALPixProjection())
        for coords in (SVector(NaN, 1.0), SVector(1.0, Inf))
            @test all(isnan, @inferred(IntermediateToNative(projection)(coords)))
            @test all(isnan, @inferred(NativeToIntermediate(projection)(coords)))
        end
    end
end

# Projections are NaN-only: far outside every domain they must return NaN in every component, never throw or overflow.
@testset "Projections never throw on invalid coordinates" begin
    projections = (GnomonicProjection(), StereographicProjection(), ZenithalEquidistant(), ZenithalEqualArea(), SlantOrthographic(), SlantOrthographic(0.1, 0.2),
        ZenithalPerspective(), ZenithalPerspective(1.5, 15.0), ZenithalPerspective(-0.5, 0.0), SlantZenithalPerspective(1.5, 15.0, 80.0), AiryProjection(), AiryProjection(75.8),
        ZenithalPolynomial((0.0, 0.95, -0.025, -0.15833, 0.00208, 0.00792, -0.00007, -0.00019)), ZenithalPolynomial((0.0, 1.0, -0.2)),
        PlateCarreeProjection(), CylindricalEqualArea(0.75), CylindricalPerspective(), CylindricalPerspective(-1.5, 0.5), MercatorProjection(), SansonFlamsteedProjection(), ParabolicProjection(),
        MollweideProjection(), PolyconicProjection(), HammerAitoffProjection(), ConicPerspective(60.0, 15.0), ConicEquidistant(45.0, 25.0), ConicEqualArea(-45.0, 10.0),
        ConicOrthomorphic(45.0, 10.0), BonneProjection(45.0), BonneProjection(0.0), TangentialSphericalCube(), COBESphericalCube(), QuadrilateralizedSphericalCube(),
        HEALPixProjection(), HEALPixProjection(5, 2), HEALPixPolarProjection())
    plane = [-1000.0, -360.0, -270.0, -180.0, -135.0, -90.0, -45.0, -1e-9, 0.0, 1e-9, 30.0, 45.0, 90.0, 135.0, 180.0, 270.0, 360.0, 1000.0]
    native_lon = [-720.0, -360.0, -180.0, -90.0, -1e-9, 0.0, 1e-9, 45.0, 90.0, 179.999999, 180.0, 360.0, 720.0]
    native_lat = [-180.0, -91.0, -90.0, -89.999999, -45.0, -1e-9, 0.0, 1e-9, 45.0, 89.999999, 90.0, 91.0, 180.0]
    problems = String[]
    for p in projections, (transform, grid) in ((IntermediateToNative(p), Iterators.product(plane, plane)), (NativeToIntermediate(p), Iterators.product(native_lon, native_lat)))
        for (a, b) in grid
            result = try
                transform(SVector(a, b))
            catch error
                push!(problems, "$(GWCS._describe_transform(transform)) $((a, b)) threw $(sprint(showerror, error))"); continue
            end
            all(isnan, result) || all(isfinite, result) || push!(problems, "$(GWCS._describe_transform(transform)) $((a, b)) returned $result")
        end
    end
    @test isempty(problems)
end

@testset "HEALPix projection parameters" begin
    # References generated with Astropy's WCSLIB-backed HEALPix models.
    references = (
        (HEALPixProjection(4, 4), SVector(30.0, -80.0), SVector(7.3954100059463315, -101.4068849910805)),
        (HEALPixProjection(4, 4), SVector(-170.0, -60.0), SVector(-172.67949192431124, -79.55771365940052)),
        (HEALPixProjection(5, 2), SVector(30.0, -80.0), SVector(34.9541310870281, -47.72478652216861)),
        (HEALPixProjection(5, 4), SVector(30.0, -80.0), SVector(34.52091799881073, -81.1255079928644))
    )
    for (projection, native, expected) in references
        forward = NativeToIntermediate(projection)
        inverse = IntermediateToNative(projection)
        @test forward(native) ≈ expected atol=3e-14
        @test inverse(forward(native)) ≈ native atol=3e-13
    end
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
    @test_throws "native reference coordinate must be a tuple or vector" celestial_transform(GnomonicProjection(), (0.0, 0.0); native_reference=0.0)
    @test_throws "native reference coordinate has length 1, expected 2" celestial_transform(GnomonicProjection(), (0.0, 0.0); native_reference=(0.0,))
    @test_throws "native reference coordinates must be real numbers" celestial_transform(GnomonicProjection(), (0.0, 0.0); native_reference=(0.0, "north"))
    @test_throws "lonpole must be a real number or nothing" celestial_transform(GnomonicProjection(), (0.0, 0.0); lonpole="north")
    @test_throws "offset must be a boolean" celestial_transform(GnomonicProjection(), (0.0, 0.0); offset=1)
    @test_throws "reference latitude must be between -90 and 90 degrees" celestial_transform(GnomonicProjection(), (0.0, 91.0))
    _test_projection_performance(transform, SVector(1.0, 2.0))
    _test_projection_performance(inverse_transform(transform), transform(SVector(1.0, 2.0)))

    offset_transform = celestial_transform(GnomonicProjection(), (10.0, 20.0); native_reference=[5.0, 30.0], offset=true)
    offset_coords = SVector(1.0, 2.0)
    offset_world = offset_transform(offset_coords)
    @test offset_transform(SVector(0.0, 0.0)) ≈ SVector(10.0, 20.0) atol=1e-12
    @test inverse_transform(offset_transform)(offset_world) ≈ offset_coords atol=1e-10
    _test_projection_performance(offset_transform, offset_coords)
    _test_projection_performance(inverse_transform(offset_transform), offset_world)

    azp_transform = celestial_transform(ZenithalPerspective(1.5, 15.0), (83.8221, -5.3911))
    wcs = WCS(Step(Frame2D(; name="intermediate"), azp_transform), Step(CelestialFrame(; name="world")))
    world = pixel_to_world(wcs, SVector(1.0, 2.0))
    @test world_to_pixel(wcs, world) ≈ SVector(1.0, 2.0) atol=1e-10
    @test @inferred(pixel_to_world(wcs, SVector(1.0, 2.0))) == world
end
