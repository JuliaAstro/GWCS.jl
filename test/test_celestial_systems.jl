using ASDF
using Logging
using GWCS: CelestialFrame, CompositeFrame, CoordinateFrame, IdentityTransform, SpectralFrame, Step, axis_names, celestial_frame, celestial_system, get_frame,
            with_celestial_system

# The angle between two positions in degrees, by Vincenty's formula, which stays accurate for tiny angles.
function _separation(l1, b1, l2, b2)
    dl, b1, b2 = deg2rad(l2 - l1), deg2rad(b1), deg2rad(b2)
    return rad2deg(atan(hypot(cos(b2) * sin(dl), cos(b1) * sin(b2) - sin(b1) * cos(b2) * cos(dl)), sin(b1) * sin(b2) + cos(b1) * cos(b2) * cos(dl)))
end

@testset "Celestial system conversion" begin
    # validation/generate_celestial_conversions_fixture.py: positions converted between each pair of systems by astropy.
    fixture = with_logger(NullLogger()) do
        ASDF.load(joinpath(@__DIR__, "fixtures", "celestial_conversions.asdf"); extensions=true)
    end
    reference(spec) = spec["system"] == "icrs" ? :ICRS : spec["system"] == "galactic" ? GWCS._fits_astropy_frame("galactic") :
        spec["system"] == "supergalactic" ? :Supergalactic : GWCS._fits_astropy_frame("fk5", "equinox" => GWCS._fits_epoch(spec["equinox"], "jyear"))
    sky(spec) = WCS(Step(CoordinateFrame{2}(; name="input"), IdentityTransform{2}()), Step(CelestialFrame(; name="world", reference_frame=reference(spec))))
    points = permutedims(hcat(fixture["points"]["lon"][], fixture["points"]["lat"][]))
    for (source, source_spec) in fixture["systems"], (target, target_spec) in fixture["systems"]
        options = target_spec["system"] == "fk5" ? (; equinox=target_spec["equinox"]) : (;)
        for keep_frame in (false, true)
            converted = with_celestial_system(sky(source_spec), Symbol(target_spec["system"]); keep_frame, options...)
            world = pixel_to_world(converted, points)
            expected = fixture["converted"][source][target]
            @test maximum(k -> _separation(world[1, k], world[2, k], expected["lon"][][k], expected["lat"][][k]), axes(points, 2)) < 1e-12
            @test all(lon -> 0 <= lon < 360, world[1, :])
            @test celestial_system(celestial_frame(converted)).system === Symbol(target_spec["system"])
            back = world_to_pixel(converted, world)
            @test maximum(k -> _separation(back[1, k], back[2, k], points[1, k], points[2, k]), axes(points, 2)) < 1e-12
        end
    end

    # mod gives 360 for a longitude just below zero, which CartesianToSpherical wraps to 0.
    @test GWCS.CartesianToSpherical()(SVector(1.0, -1e-17, 0.0))[1] === 0.0

    tan = Dict("CTYPE1" => "GLON-TAN", "CTYPE2" => "GLAT-TAN", "CRPIX1" => 50.0, "CRPIX2" => 50.0, "CRVAL1" => 120.0, "CRVAL2" => 20.0,
        "CDELT1" => -1e-3, "CDELT2" => 1e-3)
    galactic = parse_fits_wcs(tan)
    icrs = with_celestial_system(galactic, :icrs)
    kept = with_celestial_system(galactic, :icrs; keep_frame=true)
    @test available_frames(icrs) == available_frames(galactic)
    @test available_frames(kept) == (available_frames(galactic)..., "icrs")
    @test celestial_system(get_frame(kept, "world")).system === :galactic
    # The rotation is folded into the FITS NativeToCelestial, so the replaced WCS has the original's structure.
    @test typeof(icrs.forward) == typeof(galactic.forward)
    @test pixel_to_world(icrs, SVector(10.0, 70.0)) ≈ pixel_to_world(kept, SVector(10.0, 70.0)) atol=1e-12
    @test output_frame(icrs) isa CelestialFrame{2} && axis_names(output_frame(icrs)) == ("ra", "dec")
    @test world_axis_physical_types(icrs) == ("pos.eq.ra", "pos.eq.dec")
    @test celestial_system(celestial_frame(icrs)) == (system=:icrs, equinox=nothing, obstime=nothing)
    fk5 = celestial_system(celestial_frame(with_celestial_system(galactic, :fk5; equinox=1975)))
    @test fk5 == (system=:fk5, equinox=1975.0, obstime=nothing)
    @test celestial_system(celestial_frame(with_celestial_system(galactic, :fk5))).equinox == 2000.0
    # Converting back to the original system reproduces its coordinates.
    pixel = SVector(10.0, 70.0)
    @test pixel_to_world(with_celestial_system(icrs, :galactic; name="back"), pixel) ≈ pixel_to_world(galactic, pixel) atol=1e-11
    @test world_to_pixel(icrs, pixel_to_world(icrs, pixel)) ≈ pixel atol=1e-9
    @test all(isnan, pixel_to_world(icrs, SVector(NaN, 1.0)))
    allocated_forward(wcs, value) = @allocated pixel_to_world(wcs, value)
    allocated_inverse(wcs, value) = @allocated world_to_pixel(wcs, value)
    world = pixel_to_world(icrs, pixel)
    world_to_pixel(icrs, world)
    @test @inferred(pixel_to_world(icrs, pixel)) == world
    @test allocated_forward(icrs, pixel) == 0
    @test allocated_inverse(icrs, world) == 0

    # In a cube, the celestial axes keep their positions and the spectral axis passes through.
    cube = Dict("NAXIS" => 3, "CTYPE1" => "FREQ", "CTYPE2" => "GLAT-TAN", "CTYPE3" => "GLON-TAN", "CRPIX1" => 1.0, "CRPIX2" => 1.0, "CRPIX3" => 1.0,
        "CRVAL1" => 1e9, "CRVAL2" => 20.0, "CRVAL3" => 120.0, "CDELT1" => 1e6, "CDELT2" => 1e-3, "CDELT3" => 1e-3, "CUNIT1" => "Hz")
    cube_galactic = parse_fits_wcs(cube)
    cube_icrs = with_celestial_system(cube_galactic, :icrs; name="sky")
    @test output_frame(cube_icrs) isa CompositeFrame && output_frame(cube_icrs).name == "sky"
    @test celestial_frame(cube_icrs).axes_order == (3, 2)
    @test any(frame -> frame isa SpectralFrame, output_frame(cube_icrs).frames)
    voxel = SVector(4.0, 5.0, 6.0)
    original, converted = pixel_to_world(cube_galactic, voxel), pixel_to_world(cube_icrs, voxel)
    @test converted[1] == original[1]
    icrs_sky = pixel_to_world(with_celestial_system(parse_fits_wcs(merge(tan, Dict("CRVAL1" => original[3], "CRVAL2" => original[2], "CRPIX1" => 1.0, "CRPIX2" => 1.0))), :icrs), SVector(0.0, 0.0))
    @test _separation(converted[3], converted[2], icrs_sky...) < 1e-12
    @test world_to_pixel(cube_icrs, converted) ≈ voxel atol=1e-8
    # A composite output frame is not fused: the rotation step is composed after the transform into it.
    @test available_frames(cube_icrs) == (available_frames(cube_galactic)[1:(end - 1)]..., "sky")
    @test pixel_to_world(cube_icrs, voxel) ≈ pixel_to_world(with_celestial_system(cube_galactic, :icrs; keep_frame=true), voxel) atol=1e-12

    # A Python GWCS V2/V3-to-sky rotation absorbs the conversion into its Rotation3D, dropping the unmodeled properties
    # that described the old one.
    roman = parse_wcs(read(joinpath(@__DIR__, "fixtures", "roman_wcs.yaml"), String))
    roman_galactic = with_celestial_system(roman, :galactic)
    @test typeof(roman_galactic.forward) == typeof(roman.forward)
    @test pixel_to_world(roman_galactic, SVector(11.0, 22.0)) ≈ pixel_to_world(with_celestial_system(roman, :galactic; keep_frame=true), SVector(11.0, 22.0)) atol=1e-12
    bounded = GWCS.AnnotatedTransform(GWCS.NativeToCelestial(10.0, 20.0, 180.0); name="sky", pass_through=GWCS._PassThrough("" => ("tag:stsci.edu:asdf/transform/rotate3d", Dict{Any,Any}("bounds" => Dict("phi" => [0, 1])))))
    rotated = with_celestial_system(WCS(Step(CoordinateFrame{2}(; name="native"), bounded), Step(CelestialFrame(; name="world"))), :galactic)
    @test rotated.steps[1].transform isa GWCS.AnnotatedTransform && rotated.steps[1].transform.name == "sky"
    @test isempty(rotated.steps[1].transform.pass_through)

    # Written to ASDF with the new frame's astropy reference frame, and read back to the same coordinates.
    mktempdir(@__DIR__) do directory
        filename = joinpath(directory, "converted.asdf")
        fk5_wcs = with_celestial_system(galactic, :fk5; equinox=1975.0)
        write_asdf(filename, fk5_wcs)
        loaded = with_logger(NullLogger()) do
            ASDF.load(filename; extensions=true)
        end
        reparsed = parse_wcs(loaded["wcs"])
        @test pixel_to_world(reparsed, pixel) ≈ pixel_to_world(fk5_wcs, pixel) atol=1e-11
        @test celestial_system(celestial_frame(reparsed)) == (system=:fk5, equinox=1975.0, obstime=nothing)
        @test_throws "reference frame :Supergalactic has no ASDF representation" write_asdf(filename, with_celestial_system(galactic, :supergalactic))
        # A Supergalactic FITS image converted to ICRS writes once its frame is replaced, but not while it is kept.
        supergalactic = parse_fits_wcs(merge(tan, Dict("CTYPE1" => "SLON-TAN", "CTYPE2" => "SLAT-TAN")))
        write_asdf(filename, with_celestial_system(supergalactic, :icrs))
        reloaded = with_logger(NullLogger()) do
            parse_wcs(ASDF.load(filename; extensions=true)["wcs"])
        end
        @test pixel_to_world(reloaded, pixel) ≈ pixel_to_world(with_celestial_system(supergalactic, :icrs), pixel) atol=1e-11
        @test_throws "has no ASDF representation" write_asdf(filename, with_celestial_system(supergalactic, :icrs; keep_frame=true))
    end

    @test_throws "has no celestial frame" with_celestial_system(parse_fits_wcs(Dict("CTYPE1" => "FREQ", "CRPIX1" => 1.0, "CRVAL1" => 1e9, "CDELT1" => 1e6)), :icrs)
    # Other angle units are scaled to degrees around the rotation, and the new frame keeps them.
    arcsec = parse_fits_wcs(merge(tan, Dict("CUNIT1" => "arcsec", "CUNIT2" => "arcsec", "CRVAL1" => 120.0 * 3600, "CRVAL2" => 20.0 * 3600, "CDELT1" => -3.6, "CDELT2" => 3.6)); preserve_units=true)
    arcsec_icrs = with_celestial_system(arcsec, :icrs)
    @test world_axis_units(arcsec_icrs) == ("arcsec", "arcsec")
    @test pixel_to_world(arcsec_icrs, pixel) ≈ 3600 .* pixel_to_world(icrs, pixel) rtol=1e-14
    @test world_to_pixel(arcsec_icrs, pixel_to_world(arcsec_icrs, pixel)) ≈ pixel atol=1e-9
    @test_throws "has unit \"pixel\"; its axes must be in an angle unit (arcmin, arcsec, deg, hourangle, mas, rad)" with_celestial_system(WCS(Step(CoordinateFrame{2}(; name="input"), IdentityTransform{2}()), Step(CelestialFrame(; name="world", unit=("pixel", "deg")))), :icrs)
    ecliptic = parse_fits_wcs(merge(tan, Dict("CTYPE1" => "ELON-TAN", "CTYPE2" => "ELAT-TAN")))
    @test_throws "the WCS's celestial system :ecliptic cannot be converted" with_celestial_system(ecliptic, :icrs)
    @test_throws "celestial system :fk4 cannot be converted" with_celestial_system(galactic, :fk4)
    @test_throws "equinox applies only to an FK5 system" with_celestial_system(galactic, :icrs; equinox=2000.0)
    @test_throws "the WCS already has a frame named \"intermediate\"" with_celestial_system(galactic, :icrs; name="intermediate")
    @test_throws "equinox must be finite" with_celestial_system(galactic, :fk5; equinox=NaN)

    # A final remapping whose transform does not end in a sky rotation is not fused.
    remapped_sky = WCS(Step(CoordinateFrame{2}(; name="native"), GWCS.SphericalToCartesian() | GWCS.annotate_transform(GWCS.CartesianToSpherical() | GWCS.RemapAxes{2}((2, 1)); name="sky")),
        Step(CelestialFrame(; name="world", axes_order=(2, 1))))
    unit_vector = SVector(10.0, 20.0)
    @test pixel_to_world(with_celestial_system(remapped_sky, :galactic), unit_vector) ≈ pixel_to_world(with_celestial_system(remapped_sky, :galactic; keep_frame=true), unit_vector) atol=1e-12

    # Fused into a NativeToCelestial whose output is near the new pole, latitudes keep full precision.
    near_pole = parse_fits_wcs(Dict("CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CRPIX1" => 50.0, "CRPIX2" => 50.0, "CRVAL1" => 192.8594812065348,
        "CRVAL2" => 27.12825118085622, "CDELT1" => -1e-3, "CDELT2" => 1e-3))
    pole_fused, pole_kept = with_celestial_system(near_pole, :galactic), with_celestial_system(near_pole, :galactic; keep_frame=true)
    for pixel in (SVector(49.0, 49.0), SVector(49.0, 49.001), SVector(10.0, 70.0))
        @test pixel_to_world(pole_fused, pixel)[2] ≈ pixel_to_world(pole_kept, pixel)[2] atol=1e-12
        @test world_to_pixel(pole_fused, pixel_to_world(pole_fused, pixel)) ≈ pixel atol=1e-8
    end

    # Replacing a composite frame keeps the metadata of the step and of its unchanged constituents.
    notes = GWCS._PassThrough("" => ("tag:stsci.edu:gwcs/step", Dict{Any,Any}("note" => "step")),
        "frame/frames/1" => ("tag:stsci.edu:gwcs/celestial_frame", Dict{Any,Any}("note" => "sky")), "frame/frames/2" => ("tag:stsci.edu:gwcs/spectral_frame", Dict{Any,Any}("note" => "wave")))
    noted = WCS(Step(CoordinateFrame{3}(; name="pixels"), IdentityTransform{3}()),
        Step(CompositeFrame((CelestialFrame(; name="sky", axes_order=(1, 2)), SpectralFrame(; name="wave", axes_order=(3,), unit="nm")); name="world"), nothing, notes))
    @test sort!(collect(keys(last(with_celestial_system(noted, :galactic).steps).pass_through))) == ["", "frame/frames/2"]
    @test_throws "the WCS already has a frame named \"world\"" with_celestial_system(galactic, :icrs; keep_frame=true, name="world")
end
