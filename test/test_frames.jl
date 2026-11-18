using ASDF
using Logging
using GWCS: CelestialFrame, CompositeFrame, CoordinateFrame, Frame2D, IdentityTransform, PIXEL, SPATIAL, SPECTRAL, STOKES, SpectralFrame, Step,
            StokesFrame, TIME, TemporalFrame, UNDEFINED, axis_names, axis_physical_types, axis_types, axis_units, celestial_frame,
            celestial_system, fix_inputs, frame_dimension, offset_to_time, time_to_offset

@testset "Coordinate frames" begin
    generic = CoordinateFrame{3}(; name="generic", axes_names=("x", "wave", "quality"), axis_types=(PIXEL, "spectral", "custom:quality"))
    @test axis_names(generic) == ("x", "wave", "quality")
    @test axis_types(generic) == (PIXEL, SPECTRAL, "custom:quality")
    @test_throws "custom axis type must not be empty" CoordinateFrame{1}(; name="bad", axis_types=("",))
    @test_throws "axis_types has length 1, expected 2" CoordinateFrame{2}(; name="bad", axis_types=(PIXEL,))
    @test_throws "axes_names has length 1, expected 2" CoordinateFrame{2}(; name="bad", axes_names=("x",))
    @test axis_types(CoordinateFrame{1}(; name="scalar", axis_types="pixel")) == (PIXEL,)
    @test CoordinateFrame{1}(; name="scalar", axes_order=1).axes_order == (1,)
    @test CoordinateFrame{2}(; name="compatible") isa CoordinateFrame{2,Nothing}
    @test CelestialFrame(; name="compatible") isa CelestialFrame{2,Symbol}
    positional_generic = CoordinateFrame{2,Nothing}("generic", ("x", "y"), (1, 2), ("custom:x", "custom:y"), ("", ""), nothing)
    @test axis_types(positional_generic) == (UNDEFINED, UNDEFINED)
    @test CoordinateFrame("generic", ("x", "y"), (1, 2), ("custom:x", "custom:y"), ("", ""), nothing) == positional_generic
    positional_celestial = CelestialFrame{2,Symbol}("sky", ("lon", "lat"), (1, 2), ("pos.eq.ra", "pos.eq.dec"), ("deg", "deg"), :ICRS)
    @test axis_types(positional_celestial) == (SPATIAL, SPATIAL)
    @test CelestialFrame("sky", ("lon", "lat"), (1, 2), ("pos.eq.ra", "pos.eq.dec"), ("deg", "deg"), :ICRS) == positional_celestial

    pixels = Frame2D(; name="detector", unit=("pixel", "pix"))
    @test axis_types(pixels) == (SPATIAL, SPATIAL)
    @test axis_types(Frame2D(; name="focal", unit=("deg", "deg"))) == (SPATIAL, SPATIAL)

    sky = CelestialFrame(; name="sky", axes_order=(2, 5))
    @test axis_types(sky) == (SPATIAL, SPATIAL)
    @test axis_units(sky) == ("deg", "deg")

    @test axis_physical_types(SpectralFrame(; unit="GHz")) == ("em.freq",)
    @test axis_physical_types(SpectralFrame(; unit="nm")) == ("em.wl",)
    @test axis_physical_types(SpectralFrame(; unit="keV")) == ("em.energy",)
    @test axis_physical_types(SpectralFrame(; unit="km/s")) == ("spect.dopplerVeloc",)
    @test axis_physical_types(SpectralFrame(; unit="adu")) == ("custom:spectral",)
    @test axis_physical_types(SpectralFrame(; unit="nm", axis_physical_types=("em.wl;obs",))) == ("em.wl;obs",)
    @test SpectralFrame(; reference_position="barycenter").reference_position == "barycenter"
    @test_throws "unsupported spectral reference position" SpectralFrame(; reference_position="observer")

    temporal = TemporalFrame(; name="time", reference_epoch=100.0, time_scale="utc", time_format="jd", location="geocenter")
    @test axis_types(temporal) == (TIME,)
    @test offset_to_time(temporal, 86400.0) == 101.0
    @test time_to_offset(temporal, 101.0) == 86400.0
    @test temporal.time_scale == "utc"
    @test temporal.time_format == "jd"
    @test temporal.location == "geocenter"
    @test_throws "has no reference epoch" offset_to_time(TemporalFrame(), 1.0)
    @test_throws "cannot convert temporal unit fortnight to days" offset_to_time(TemporalFrame(; reference_epoch=100.0, time_format="jd", unit="fortnight"), 1.0)
    @test_throws "numeric temporal format jyear requires an optional time integration" offset_to_time(TemporalFrame(; reference_epoch=2000.0, time_format="jyear"), 86400.0)

    stokes = StokesFrame()
    @test axis_types(stokes) == (STOKES,)
    @test axis_physical_types(stokes) == ("phys.polarization.stokes",)
    @test axis_units(stokes) == ("",)

    spectral = SpectralFrame(; name="spectral", axes_order=(1,), unit="nm")
    temporal = TemporalFrame(; name="temporal", axes_order=(4,), reference_epoch=0.0)
    stokes = StokesFrame(; name="polarization", axes_order=(3,))
    composite = CompositeFrame((sky, spectral, temporal, stokes); name="world")
    @test frame_dimension(composite) == 5
    @test composite.axes_order == (2, 5, 1, 4, 3)
    # The accessors follow the transform order: axes_order puts the spectral axis first and the latitude last.
    @test axis_names(composite) == ("spectral", "lon", "stokes", "time", "lat")
    @test axis_types(composite) == (SPECTRAL, SPATIAL, STOKES, TIME, SPATIAL)
    @test axis_units(composite) == ("nm", "deg", "", "s", "deg")
    @test composite.axes_names == ("lon", "lat", "spectral", "time", "stokes")
    @test composite.frames == (sky, spectral, temporal, stokes)

    @test_throws "requires at least one constituent" CompositeFrame(())
    repeated_names = CompositeFrame((SpectralFrame(; name="spectral", axes_order=1), SpectralFrame(; name="spectral", axes_order=2)))
    @test length(repeated_names.frames) == 2
    @test_throws "axes_order entries must be unique" CompositeFrame((CelestialFrame(; name="a", axes_order=(1, 2)), SpectralFrame(; name="b", axes_order=(2,))))
    @test_throws "must be a permutation of 1:3" CompositeFrame((CelestialFrame(; name="a", axes_order=(1, 2)), SpectralFrame(; name="b", axes_order=(4,))))

    nested = CompositeFrame((CompositeFrame((CelestialFrame(; name="nested_sky", axes_order=(1, 2)),); name="inner"), SpectralFrame(; name="nested_wave", axes_order=3)); name="outer")
    @test length(nested.frames) == 2
    @test all(frame -> !(frame isa CompositeFrame), nested.frames)
    subset_composite(frame) = GWCS._subset_frame(frame, (1, 3))
    @test @inferred(subset_composite(CompositeFrame((CelestialFrame(; name="inferred_sky", axes_order=(1, 2)), SpectralFrame(; name="inferred_wave", axes_order=3))))) isa CompositeFrame{2}
end

@testset "Composite frame WCS" begin
    input = CoordinateFrame{5}(; name="input", axes_names=("a", "b", "c", "d", "e"), axis_types=ntuple(_ -> PIXEL, 5))
    sky = CelestialFrame(; name="sky", axes_order=(2, 5))
    spectral = SpectralFrame(; name="spectral", axes_order=(1,), axes_names=("wave",), unit="nm")
    temporal = TemporalFrame(; name="temporal", axes_order=(4,), reference_epoch=100.0)
    stokes = StokesFrame(; name="polarization", axes_order=(3,))
    world = CompositeFrame((sky, spectral, temporal, stokes); name="world")
    wcs = WCS(Step(input, IdentityTransform{5}()), Step(world))

    # World coordinates are the transform outputs; axes_order only labels them.
    pixels = SVector(1.0, 2.0, 3.0, 4.0, 5.0)
    coordinates = pixels
    @test pixel_to_world(wcs, pixels) == coordinates
    @test world_to_pixel(wcs, coordinates) == pixels
    @test pixel_axis_names(wcs) == ("a", "b", "c", "d", "e")
    @test world_axis_names(wcs) == ("wave", "lon", "stokes", "time", "lat")
    @test world_axis_physical_types(wcs) == ("em.wl", "pos.eq.ra", "phys.polarization.stokes", "time", "pos.eq.dec")
    @test world_axis_units(wcs) == ("nm", "deg", "", "s", "deg")
    @test world_axis_types(wcs) == (SPECTRAL, SPATIAL, STOKES, TIME, SPATIAL)
    @test axis_correlation_matrix(wcs) == [i == j for i in 1:5, j in 1:5]
    @test @inferred(pixel_to_world(wcs, pixels)) == coordinates
    allocated_composite_forward(wcs, coords) = @allocated pixel_to_world(wcs, coords)
    allocated_composite_backward(wcs, coords) = @allocated world_to_pixel(wcs, coords)
    pixel_to_world(wcs, pixels)
    world_to_pixel(wcs, coordinates)
    @test allocated_composite_forward(wcs, pixels) == 0
    @test allocated_composite_backward(wcs, coordinates) == 0

    sliced = fix_inputs(WCS(Step(world, IdentityTransform{5}()), Step(input)), "wave" => 1.0)
    @test input_frame(sliced) isa CompositeFrame{4}
    @test input_frame(sliced).axes_order == (1, 4, 3, 2)
    @test axis_names(input_frame(sliced)) == ("lon", "stokes", "time", "lat")
    @test length(input_frame(sliced).frames) == 3
    @test pixel_to_world(sliced, 2.0, 3.0, 4.0, 5.0) == SVector(1.0, 2.0, 3.0, 4.0, 5.0)

    inner = CompositeFrame((CelestialFrame(; name="nested_sky", axes_order=(1, 2)),); name="inner")
    nested_input = CompositeFrame((inner, SpectralFrame(; name="nested_wave", axes_order=3, axes_names="nested_wave")); name="nested")
    nested_wcs = WCS(Step(nested_input, IdentityTransform{3}()), Step(CoordinateFrame{3}(; name="output")))
    nested_slice = fix_inputs(nested_wcs, "nested_wave" => 3.0)
    @test input_frame(nested_slice) isa CompositeFrame{2}
    @test length(input_frame(nested_slice).frames) == 1
    @test pixel_to_world(nested_slice, 1.0, 2.0) == SVector(1.0, 2.0, 3.0)

end

@testset "Celestial systems" begin
    tan = Dict("CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CRPIX1" => 1.0, "CRPIX2" => 1.0, "CRVAL1" => 10.0, "CRVAL2" => 20.0, "CDELT1" => 1e-3, "CDELT2" => 1e-3)
    fits_system(extra) = celestial_system(celestial_frame(parse_fits_wcs(merge(tan, Dict(extra)))))
    # FITS RADESYS and EQUINOX defaults are applied when the header is read.
    @test fits_system(()) == (system=:icrs, equinox=nothing, obstime=nothing)
    @test fits_system(("EQUINOX" => 1950.0,)) == (system=:fk4, equinox=1950.0, obstime=1950.0)
    @test fits_system(("RADESYS" => "FK5",)) == (system=:fk5, equinox=2000.0, obstime=nothing)
    @test fits_system(("RADESYS" => "FK4-NO-E", "EQUINOX" => 1975.0)) == (system=:fk4noeterms, equinox=1975.0, obstime=1975.0)
    @test fits_system(("RADESYS" => "GAPPT", "EQUINOX" => 2000.0)) == (system=:gappt, equinox=nothing, obstime=nothing)
    @test fits_system(("CTYPE1" => "GLON-TAN", "CTYPE2" => "GLAT-TAN")) == (system=:galactic, equinox=nothing, obstime=nothing)
    # An ecliptic EQUINOX is Besselian with an FK4 RADESYS, which EQUINOX = 1950 implies, and is reported as a Julian epoch.
    @test fits_system(("CTYPE1" => "ELON-TAN", "CTYPE2" => "ELAT-TAN", "EQUINOX" => 1950.0)).equinox ≈ 1949.9997904423 atol=1e-9
    @test fits_system(("CTYPE1" => "ELON-TAN", "CTYPE2" => "ELAT-TAN", "EQUINOX" => 2000.0)) == (system=:ecliptic, equinox=2000.0, obstime=nothing)
    @test fits_system(("CTYPE1" => "HLON-TAN", "CTYPE2" => "HLAT-TAN")) == (system=:helioecliptic, equinox=nothing, obstime=nothing)
    @test fits_system(("CTYPE1" => "SLON-TAN", "CTYPE2" => "SLAT-TAN")) == (system=:supergalactic, equinox=nothing, obstime=nothing)
    generic = celestial_frame(parse_fits_wcs(merge(tan, Dict("CTYPE1" => "XLAT-TAN", "CTYPE2" => "XLON-TAN"))))
    @test celestial_system(generic) == (system=:unknown, equinox=nothing, obstime=nothing)
    @test generic.axes_order == (2, 1)

    # In a cube, axes_order gives the celestial axes' positions in the world coordinates.
    cube = Dict("NAXIS" => 3, "CTYPE1" => "FREQ", "CTYPE2" => "RA---TAN", "CTYPE3" => "DEC--TAN", "CRPIX1" => 1.0, "CRPIX2" => 1.0, "CRPIX3" => 1.0,
        "CRVAL1" => 1e9, "CRVAL2" => 10.0, "CRVAL3" => 20.0, "CDELT1" => 1e6, "CDELT2" => 1e-3, "CDELT3" => 1e-3, "CUNIT1" => "Hz")
    @test celestial_frame(parse_fits_wcs(cube)).axes_order == (2, 3)
    @test celestial_frame(parse_fits_wcs(Dict("CTYPE1" => "FREQ", "CRPIX1" => 1.0, "CRVAL1" => 1e9, "CDELT1" => 1e6))) === nothing
    @test celestial_frame(SpectralFrame(; unit="nm")) === nothing
    twice = CompositeFrame((CelestialFrame(; name="a", axes_order=(1, 2)), CelestialFrame(; name="b", axes_order=(3, 4))); name="both")
    @test_throws "frame both has 2 celestial frames" celestial_frame(twice)

    @test celestial_system(CelestialFrame(; name="sky")) == (system=:icrs, equinox=nothing, obstime=nothing)
    @test celestial_system(CelestialFrame(; name="sky", reference_frame=:FK5)) == (system=:fk5, equinox=2000.0, obstime=nothing)
    @test celestial_system(CelestialFrame(; name="sky", reference_frame=nothing)) == (system=:unknown, equinox=nothing, obstime=nothing)
    time(value, format, scale...; version="1.4.0") = GWCS.TaggedMapping("tag:stsci.edu:asdf/time/time-$version", Dict{Any,Any}("value" => value, "format" => format, (isempty(scale) ? () : ("scale" => only(scale),))...))
    frame(name, attributes...; version="1.2.0") = CelestialFrame(; name="sky", reference_frame=GWCS.TaggedMapping("tag:astropy.org:astropy/coordinates/frames/$name-$version", Dict{Any,Any}("frame_attributes" => Dict{Any,Any}(attributes...))))
    @test celestial_system(frame("fk5", "equinox" => time(2451545.0, "jd", "tt"))).equinox == 2000.0
    @test celestial_system(frame("fk5", "equinox" => time(51544.5, "mjd", "tdb"))).equinox == 2000.0
    @test celestial_system(frame("fk5")).equinox == 2000.0
    # An FK4 obstime defaults to the equinox.
    @test celestial_system(frame("fk4", "equinox" => time(1975.0, "byear", "tt"))) == (system=:fk4, equinox=1975.0, obstime=1975.0)
    @test_throws "time scale \"utc\" is not supported" celestial_system(frame("fk5", "equinox" => time(2451545.0, "jd")))
    @test_throws "with format \"unix\" is not supported" celestial_system(frame("fk5", "equinox" => time(0.0, "unix", "tt")))
    @test_throws "time \"J1975.0\" does not match its format \"iso\"" celestial_system(frame("fk5", "equinox" => time("J1975.0", "iso", "tt")))
    @test celestial_system(frame("fk5", "equinox" => time("J1975.0", "jyear_str"))).equinox == 1975.0
    @test_throws "unsupported metadata tag version tag:stsci.edu:asdf/time/time-1.5.0" celestial_system(frame("fk5", "equinox" => time(2000.0, "jyear"; version="1.5.0")))
    @test_throws "an epoch must be a tagged ASDF time" celestial_system(frame("fk5", "equinox" => 2000.0))
    @test_throws "unsupported metadata tag version tag:astropy.org:astropy/coordinates/frames/fk5-1.3.0" celestial_system(frame("fk5"; version="1.3.0"))
    # Frames whose attributes are not read report their system at any version.
    @test celestial_system(frame("galactic"; version="9.0.0")).system === :galactic

    # validation/generate_celestial_systems_fixture.py: Python GWCS frames, with the system and epochs astropy reports.
    fixture = with_logger(NullLogger()) do
        ASDF.load(joinpath(@__DIR__, "fixtures", "celestial_systems.asdf"); extensions=true)
    end
    for (name, reference) in fixture["references"]
        system = celestial_system(celestial_frame(parse_wcs(fixture["wcs"][name])))
        @test system.system === Symbol(reference["system"])
        @test system.equinox === nothing ? !haskey(reference, "equinox") : isapprox(system.equinox, reference["equinox"]; atol=1e-9)
        @test system.obstime === nothing ? !haskey(reference, "obstime") : isapprox(system.obstime, reference["obstime"]; atol=1e-9)
    end
    # asdf-astropy writes an obstime given as a date in UTC, which needs leap seconds.
    utc = celestial_frame(parse_wcs(fixture["wcs"]["fk4noeterms_utc"]))
    @test_throws "only J and B epoch strings are read, not dates at frame world/reference_frame/frame_attributes/obstime" celestial_system(utc)
end
