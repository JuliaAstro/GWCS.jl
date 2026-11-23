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
    function reference(spec)
        spec["system"] == "icrs" && return :ICRS
        spec["system"] == "galactic" && return GWCS._fits_astropy_frame("galactic")
        spec["system"] == "supergalactic" && return :Supergalactic
        spec["system"] == "fk5" && return GWCS._fits_astropy_frame("fk5", "equinox" => GWCS._fits_epoch(spec["equinox"], "jyear"))
        # A FITS ecliptic, as parse_fits_wcs keeps it.
        spec["system"] == "ecliptic_fk5" && return (system=("ELON", "ELAT"), radesys="FK5", equinox=spec["equinox"], obstime=nothing)
        spec["system"] == "ecliptic_fk4" && return (system=("ELON", "ELAT"), radesys="FK4", equinox=spec["equinox"], obstime=spec["obstime"])
        spec["system"] == "ecliptic_icrs" && return (system=("ELON", "ELAT"), radesys="ICRS", equinox=nothing, obstime=nothing)
        return GWCS._fits_astropy_frame(spec["system"], "equinox" => GWCS._fits_epoch(spec["equinox"], "byear"), "obstime" => GWCS._fits_epoch(spec["obstime"], "byear"))
    end
    sky(spec) = WCS(Step(CoordinateFrame{2}(; name="input"), IdentityTransform{2}()), Step(CelestialFrame(; name="world", reference_frame=reference(spec))))
    points = permutedims(hcat(fixture["points"]["lon"][], fixture["points"]["lat"][]))
    fk4(spec) = spec["system"] in ("fk4", "fk4noeterms", "ecliptic_fk4")
    galactic(spec) = spec["system"] in ("galactic", "supergalactic")
    # Every conversion from or to FK4 passes through FK5 J2000. Murray's rotation, made orthogonal, differs from
    # astropy's by up to 0.1 µas; FK5 to FK4 is the exact inverse of FK4 to FK5, where astropy evaluates Newcomb's
    # precession in reverse (6 µas to an FK4 equinox of B1975, 36 µas to B2010).
    function tolerance(source, target)
        fk4(target) && target["equinox"] != 1950 && return 2e-8
        fk4(source) || fk4(target) ? 1e-10 : 1e-12
    end
    # astropy's own transform_to converts FK4 to Galactic coordinates through the B1950 pole and between FK4 systems by
    # precession alone, neither with Murray's obstime correction, so it is compared only where the obstimes cannot
    # matter: an FK4 obstime of B1950 against Galactic coordinates, and equal obstimes between FK4 systems.
    function comparable(source, target)
        fk4(source) && fk4(target) && return source["obstime"] == target["obstime"]
        fk4(source) && galactic(target) && return source["obstime"] == 1950
        fk4(target) && galactic(source) && return target["obstime"] == 1950
        return true
    end
    for (source, source_spec) in fixture["systems"], (target, target_spec) in fixture["systems"]
        options = target_spec["system"] in ("fk5", "ecliptic_fk5") ? (; equinox=target_spec["equinox"]) :
            target_spec["system"] in ("fk4", "fk4noeterms", "ecliptic_fk4") ? (; equinox=target_spec["equinox"], obstime=target_spec["obstime"]) : (;)
        for keep_frame in (false, true)
            # An ICRS ecliptic warns on each conversion; the warning is tested below.
            converted = with_logger(NullLogger()) do
                with_celestial_system(sky(source_spec), Symbol(target_spec["system"]); keep_frame, options...)
            end
            world = pixel_to_world(converted, points)
            expected = fixture["converted"][source][target]
            @test maximum(k -> _separation(world[1, k], world[2, k], expected["lon"][][k], expected["lat"][][k]), axes(points, 2)) < tolerance(source_spec, target_spec)
            if haskey(fixture["default"][source], target) && comparable(source_spec, target_spec)
                # astropy's transform_to from FK4 to FK5 passes through the FK5 equinox, up to about 60 µas away.
                default = fixture["default"][source][target]
                @test maximum(k -> _separation(world[1, k], world[2, k], default["lon"][][k], default["lat"][][k]), axes(points, 2)) < 2e-8
            end
            @test all(lon -> 0 <= lon < 360, world[1, :])
            @test celestial_system(celestial_frame(converted)).system === Symbol(target_spec["system"])
            # The E-terms' addition inverts their removal to third order in their 1.7e-6 rad size, so round trips are exact.
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
    gappt_ecliptic = parse_fits_wcs(merge(tan, Dict("CTYPE1" => "ELON-TAN", "CTYPE2" => "ELAT-TAN", "RADESYS" => "GAPPT")))
    @test_throws "the WCS's celestial system :ecliptic_gappt cannot be converted" with_celestial_system(gappt_ecliptic, :icrs)
    @test_throws "the WCS's celestial system :helioecliptic cannot be converted" with_celestial_system(parse_fits_wcs(merge(tan, Dict("CTYPE1" => "HLON-TAN", "CTYPE2" => "HLAT-TAN"))), :icrs)
    @test_throws "celestial system :gcrs cannot be converted" with_celestial_system(galactic, :gcrs)
    @test_throws "obstime applies only to FK4, FK4-NO-E, and FK4 ecliptic systems" with_celestial_system(galactic, :fk5; obstime=2000.0)
    @test_throws "obstime applies only to FK4, FK4-NO-E, and FK4 ecliptic systems" with_celestial_system(galactic, :ecliptic_fk5; obstime=1950.0)
    @test_throws "equinox applies only to FK5 and FK4 systems" with_celestial_system(galactic, :icrs; equinox=2000.0)
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

    # FK4-NO-E is a rotation, so it fuses and writes; FK4's E-terms are a separate step with no ASDF form.
    fk4_header(radesys) = parse_fits_wcs(Dict("CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CRPIX1" => 50.0, "CRPIX2" => 50.0, "CRVAL1" => 150.0,
        "CRVAL2" => 30.0, "CDELT1" => -1e-3, "CDELT2" => 1e-3, "RADESYS" => radesys, "EQUINOX" => 1950.0))
    noe_icrs = with_celestial_system(fk4_header("FK4-NO-E"), :icrs)
    @test typeof(noe_icrs.forward) == typeof(fk4_header("FK4-NO-E").forward)
    fk4_icrs = with_celestial_system(fk4_header("FK4"), :icrs)
    @test pixel_to_world(fk4_icrs, pixel) ≈ pixel_to_world(with_celestial_system(fk4_header("FK4"), :icrs; keep_frame=true), pixel) atol=1e-12
    @test world_to_pixel(fk4_icrs, pixel_to_world(fk4_icrs, pixel)) ≈ pixel atol=1e-7
    fk4_world = pixel_to_world(fk4_icrs, pixel)
    world_to_pixel(fk4_icrs, fk4_world)
    @test allocated_forward(fk4_icrs, pixel) == 0
    @test allocated_inverse(fk4_icrs, fk4_world) == 0
    @test celestial_system(celestial_frame(with_celestial_system(galactic, :fk4))) == (system=:fk4, equinox=1950.0, obstime=1950.0)
    @test celestial_system(celestial_frame(with_celestial_system(galactic, :fk4noeterms; equinox=1975, obstime=2000))) == (system=:fk4noeterms, equinox=1975.0, obstime=2000.0)
    @test world_axis_physical_types(with_celestial_system(galactic, :fk4)) == ("pos.eq.ra", "pos.eq.dec")
    mktempdir(@__DIR__) do directory
        filename = joinpath(directory, "fk4.asdf")
        write_asdf(filename, noe_icrs)
        @test_throws "no metadata serializer is registered" write_asdf(filename, fk4_icrs)
    end

    # FITS Paper II's example 2 (Calabretta & Greisen 2002, sect. 7.3.2, Tables 7 and 8): a Galactic conic image whose
    # alternate description is in FK5 J2000 ecliptic coordinates. Its tabulated values are rounded to 1e-7 degrees.
    paper2 = Dict{String,Any}("NAXIS" => 2, "NAXIS1" => 2048, "NAXIS2" => 2048, "CTYPE1" => "GLON-COE", "CTYPE2" => "GLAT-COE", "PV2_1" => -25.0,
        "CRVAL1" => 90.0, "CRVAL2" => -25.0, "CTYPE1A" => "ELON-COE", "CTYPE2A" => "ELAT-COE", "PV2_1A" => -25.0, "CRVAL1A" => -7.0300934,
        "CRVAL2A" => 34.8474143, "LONPOLEA" => 6.3839706, "LATPOLEA" => 29.8114400, "RADESYSA" => "FK5")
    for suffix in ("", "A")
        merge!(paper2, Dict("CRPIX1$suffix" => 1024.5, "CRPIX2$suffix" => -1023.5, "PC1_1$suffix" => 1.0, "PC1_2$suffix" => -0.004, "PC2_1$suffix" => -0.002,
            "PC2_2$suffix" => 1.0, "CDELT1$suffix" => -0.005, "CDELT2$suffix" => 0.005))
    end
    paper2_pixel = SVector(1957.2, 775.4) .- 1
    galactic_image, ecliptic_image = parse_fits_wcs(paper2), parse_fits_wcs(paper2; key='A')
    @test celestial_system(celestial_frame(ecliptic_image)) == (system=:ecliptic_fk5, equinox=2000.0, obstime=nothing)
    @test pixel_to_world(galactic_image, paper2_pixel) ≈ SVector(85.2439814, -15.8973800) atol=1e-7
    @test _separation(pixel_to_world(ecliptic_image, paper2_pixel)..., -14.7066741 + 360, 43.0457292) < 2e-7
    converted_image = with_celestial_system(galactic_image, :ecliptic_fk5)
    @test world_axis_physical_types(converted_image) == ("pos.ecliptic.lon", "pos.ecliptic.lat")
    @test celestial_system(celestial_frame(converted_image)) == (system=:ecliptic_fk5, equinox=2000.0, obstime=nothing)
    # The conversion reproduces the paper's ecliptic coordinates to within about a milliarcsecond.
    @test _separation(pixel_to_world(converted_image, paper2_pixel)..., -14.7066741 + 360, 43.0457292) < 1.5 / 3.6e6
    @test _separation(pixel_to_world(converted_image, SVector(1024.5, -1023.5) .- 1)..., -7.0300934 + 360, 34.8474143) < 1.5 / 3.6e6
    mktempdir(@__DIR__) do directory
        @test_throws "has no ASDF representation" write_asdf(joinpath(directory, "ecliptic.asdf"), converted_image)
    end
    # ICRS defines no ecliptic, so converting from or to astropy's J2000 mean ecliptic warns.
    icrs_ecliptic = parse_fits_wcs(merge(tan, Dict("CTYPE1" => "ELON-TAN", "CTYPE2" => "ELAT-TAN")))
    @test celestial_system(celestial_frame(icrs_ecliptic)).system === :ecliptic_icrs
    @test_logs (:warn, r"ICRS defines no ecliptic") with_celestial_system(icrs_ecliptic, :icrs)
    @test_logs (:warn, r"ICRS defines no ecliptic") with_celestial_system(galactic, :ecliptic_icrs)
    @test_logs with_celestial_system(galactic, :ecliptic_fk4; equinox=1975.0)
    @test celestial_system(celestial_frame(with_celestial_system(galactic, :ecliptic_fk4; equinox=1975.0))) == (system=:ecliptic_fk4, equinox=1975.0, obstime=1975.0)
    @test celestial_system(celestial_frame(with_celestial_system(galactic, :ecliptic_fk4; obstime=1985.0))) == (system=:ecliptic_fk4, equinox=1950.0, obstime=1985.0)
    # Murray's obstime correction enters every FK4 conversion, so the route does not matter: FK4-NO-E at obstime B1990
    # reaches ICRS directly and through Galactic coordinates the same way.
    noe1990 = WCS(Step(CoordinateFrame{2}(; name="input"), IdentityTransform{2}()), Step(CelestialFrame(; name="world",
        reference_frame=GWCS._fits_astropy_frame("fk4noeterms", "equinox" => GWCS._fits_epoch(1950.0, "byear"), "obstime" => GWCS._fits_epoch(1990.0, "byear")))))
    @test pixel_to_world(with_celestial_system(noe1990, :icrs), SVector(120.0, 30.0)) ≈ pixel_to_world(with_celestial_system(with_celestial_system(noe1990, :galactic), :icrs), SVector(120.0, 30.0)) atol=1e-12
    # An ecliptic named by Symbol takes its system's default equinox.
    symbolic(system) = WCS(Step(CoordinateFrame{2}(; name="input"), IdentityTransform{2}()), Step(CelestialFrame(; name="world", reference_frame=system)))
    @test celestial_system(celestial_frame(symbolic(:ecliptic_fk5))) == (system=:ecliptic_fk5, equinox=2000.0, obstime=nothing)
    @test pixel_to_world(with_celestial_system(symbolic(:ecliptic_fk4), :icrs), pixel) ≈ pixel_to_world(with_celestial_system(symbolic(:ecliptic_fk4), :icrs; keep_frame=true), pixel) atol=1e-12

    # Replacing a composite frame keeps the metadata of the step and of its unchanged constituents.
    notes = GWCS._PassThrough("" => ("tag:stsci.edu:gwcs/step", Dict{Any,Any}("note" => "step")),
        "frame/frames/1" => ("tag:stsci.edu:gwcs/celestial_frame", Dict{Any,Any}("note" => "sky")), "frame/frames/2" => ("tag:stsci.edu:gwcs/spectral_frame", Dict{Any,Any}("note" => "wave")))
    noted = WCS(Step(CoordinateFrame{3}(; name="pixels"), IdentityTransform{3}()),
        Step(CompositeFrame((CelestialFrame(; name="sky", axes_order=(1, 2)), SpectralFrame(; name="wave", axes_order=(3,), unit="nm")); name="world"), nothing, notes))
    @test sort!(collect(keys(last(with_celestial_system(noted, :galactic).steps).pass_through))) == ["", "frame/frames/2"]
    @test_throws "the WCS already has a frame named \"world\"" with_celestial_system(galactic, :icrs; keep_frame=true, name="world")
end
