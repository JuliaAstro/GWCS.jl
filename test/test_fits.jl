import FITSFiles
import FITSIO

@testset "Native FITS WCS parsing" begin
    linear_header = Dict(
        "NAXIS" => 2, "NAXIS1" => 100, "NAXIS2" => 200,
        "CTYPE1" => "X", "CTYPE2" => "Y", "CNAME1" => "detector x",
        "CRPIX1" => 10.0, "CRPIX2" => 20.0,
        "CRVAL1" => 100.0, "CRVAL2" => 200.0,
        "CDELT1" => 2.0, "CDELT2" => 3.0,
    )
    linear = parse_fits_wcs(linear_header)
    @test linear isa WCS
    @test available_frames(linear) == ("detector", "intermediate", "world")
    @test linear.pixel_shape == (100, 200)
    @test linear.bounding_box === nothing
    @test world_axis_names(linear) == ("detector x", "Y")
    @test pixel_to_world(linear, 10.0, 20.0) ≈ SVector(100.0, 200.0)
    @test pixel_to_world(linear, 11.0, 22.0) ≈ SVector(102.0, 206.0)
    @test world_to_pixel(linear, 102.0, 206.0) ≈ SVector(11.0, 22.0)

    pc_header = merge(linear_header, Dict("PC1_1" => 1.0, "PC1_2" => 0.5, "PC2_1" => -0.25, "PC2_2" => 1.0))
    @test pixel_to_world(parse_fits_wcs(pc_header), 11.0, 22.0) ≈ SVector(104.0, 205.25)
    # Like WCSLIB and astropy, PCi_j and CDELTi take precedence over CDi_j in a nonconforming header.
    mixed_matrix = @test_logs (:warn, r"both PCi_j and CDi_j") parse_fits_wcs(merge(pc_header, Dict("CD1_1" => 7.0, "CD2_2" => 9.0)))
    @test pixel_to_world(mixed_matrix, 11.0, 22.0) ≈ SVector(104.0, 205.25)

    # Like WCSLIB, CROTA applies only to a celestial longitude/latitude pair.
    crota_header = merge(linear_header, Dict("CDELT1" => 1.0, "CDELT2" => 1.0, "CROTA2" => 90.0))
    @test pixel_to_world(parse_fits_wcs(crota_header), 11.0, 20.0) == SVector(101.0, 200.0)
    @test_throws "singular transformation matrix because CDELT1 is zero" parse_fits_wcs(Dict("NAXIS" => 2, "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CDELT1" => 0.0, "CROTA2" => 10.0))
    reordered_crota_header = Dict(
        "NAXIS" => 3, "CTYPE1" => "FREQ", "CTYPE2" => "DEC--CAR", "CTYPE3" => "RA---CAR",
        "CRPIX1" => 1.0, "CRPIX2" => 1.0, "CRPIX3" => 1.0, "CRVAL1" => 10.0, "CRVAL2" => 20.0, "CRVAL3" => 30.0,
        "CDELT1" => 2.0, "CDELT2" => 3.0, "CDELT3" => -4.0, "CROTA2" => 90.0,
    )
    reordered_crota = parse_fits_wcs(reordered_crota_header)
    @test pixel_to_world(reordered_crota, 2.0, 1.0, 1.0) ≈ SVector(12.0, 20.0, 30.0) atol=1e-13
    @test pixel_to_world(reordered_crota, 1.0, 2.0, 1.0) ≈ SVector(10.0, 19.97142297, 26.80785288) atol=5e-9
    swapped_crota_header = merge(reordered_crota_header, Dict("CTYPE2" => "RA---CAR", "CTYPE3" => "DEC--CAR", "CRVAL2" => 30.0, "CRVAL3" => 20.0, "CROTA2" => 0.0, "CROTA3" => 90.0))
    @test pixel_to_world(parse_fits_wcs(swapped_crota_header), 1.0, 2.0, 1.0) ≈ SVector(10.0, 30.0, 23.0) atol=5e-13

    tan_header = Dict(
        "NAXIS" => 2, "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN",
        "CRPIX1" => 512.0, "CRPIX2" => 512.0,
        "CRVAL1" => 83.8221, "CRVAL2" => -5.3911,
        "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4,
    )
    tan = parse_fits_wcs(tan_header)
    reference = SVector(83.8221, -5.3911)
    @test pixel_to_world(tan, 512.0, 512.0) ≈ reference atol=1e-12
    @test world_to_pixel(tan, reference) ≈ SVector(512.0, 512.0) atol=2e-11
    sample_world = pixel_to_world(tan, 400.0, 300.0)
    @test world_to_pixel(tan, sample_world) ≈ SVector(400.0, 300.0) atol=1e-9
    @test eltype(pixel_to_world(tan, SVector(400.0f0, 300.0f0))) === Float32
    pixel_to_world(tan, SVector(400.0, 300.0))
    world_to_pixel(tan, sample_world)
    allocated_fits_forward(wcs, pixel) = @allocated pixel_to_world(wcs, pixel)
    allocated_fits_inverse(wcs, world) = @allocated world_to_pixel(wcs, world)
    @test allocated_fits_forward(tan, SVector(400.0, 300.0)) == 0
    @test allocated_fits_inverse(tan, sample_world) == 0
    # Non-finite coordinates propagate as NaN rather than raising projection domain errors.
    @test all(isnan, pixel_to_world(tan, NaN, 300.0))
    @test all(isnan, world_to_pixel(tan, NaN, -5.0))
    @test all(isnan, pixel_to_world(tan, [NaN 400.0; 300.0 300.0])[:, 1])

    interleaved_header = Dict(
        "NAXIS" => 3, "CTYPE1" => "FREQ", "CTYPE2" => "DEC--CAR", "CTYPE3" => "RA---CAR",
        "CUNIT1" => "GHz", "CUNIT2" => "deg", "CUNIT3" => "deg",
        "CRPIX1" => 1.0, "CRPIX2" => 2.0, "CRPIX3" => 3.0,
        "CRVAL1" => 1.5, "CRVAL2" => -20.0, "CRVAL3" => 120.0,
        "CDELT1" => 0.1, "CDELT2" => 0.2, "CDELT3" => -0.2,
    )
    interleaved = parse_fits_wcs(interleaved_header)
    @test pixel_to_world(interleaved, 1.0, 2.0, 3.0) ≈ SVector(1.5e9, -20.0, 120.0) atol=1e-8
    @test world_axis_types(interleaved) == (SPECTRAL, SPATIAL, SPATIAL)
    @test world_axis_units(interleaved) == ("Hz", "deg", "deg")
    @test world_to_pixel(interleaved, 1.5e9, -20.0, 120.0) ≈ SVector(1.0, 2.0, 3.0) atol=1e-12

    arcsec_header = merge(tan_header, Dict(
        "CRVAL1" => tan_header["CRVAL1"] * 3600, "CRVAL2" => tan_header["CRVAL2"] * 3600,
        "CDELT1" => tan_header["CDELT1"] * 3600, "CDELT2" => tan_header["CDELT2"] * 3600,
        "CUNIT1" => "arcsec", "CUNIT2" => "arcsec",
    ))
    canonical = parse_fits_wcs(arcsec_header)
    preserved = parse_fits_wcs(arcsec_header; preserve_units=true)
    @test pixel_to_world(canonical, 512.0, 512.0) ≈ reference atol=1e-12
    @test pixel_to_world(preserved, 512.0, 512.0) ≈ reference .* 3600 atol=1e-8
    @test world_axis_units(canonical) == ("deg", "deg")
    @test world_axis_units(preserved) == ("arcsec", "arcsec")
    @test output_frame(tan).reference_frame.radesys == "ICRS"
    fk4 = parse_fits_wcs(merge(tan_header, Dict("EQUINOX" => 1950.0)))
    @test output_frame(fk4).reference_frame.radesys == "FK4"
    @test output_frame(fk4).reference_frame.equinox == 1950.0
    @test output_frame(parse_fits_wcs(merge(tan_header, Dict("RADESYS" => "FK5")))).reference_frame.equinox == 2000.0
    @test output_frame(parse_fits_wcs(merge(tan_header, Dict("RADESYS" => "FK4-NO-E")))).reference_frame.equinox == 1950.0
    icrs = output_frame(parse_fits_wcs(merge(tan_header, Dict("RADESYS" => "icrs", "EQUINOX" => 1950.0)))).reference_frame
    @test icrs.radesys == "ICRS"
    @test icrs.equinox === nothing
    gappt = output_frame(parse_fits_wcs(merge(tan_header, Dict("RADESYS" => "GAPPT", "EQUINOX" => 2000.0)))).reference_frame
    @test gappt.equinox === nothing
    ecliptic_header = merge(tan_header, Dict("CTYPE1" => "ELON-AIT", "CTYPE2" => "ELAT-AIT", "EQUINOX" => 2000.0))
    ecliptic = output_frame(parse_fits_wcs(ecliptic_header)).reference_frame
    @test ecliptic.radesys == "FK5"
    @test ecliptic.equinox == 2000.0
    generic_header = merge(tan_header, Dict("CTYPE1" => "XLAT-BON", "CTYPE2" => "XLON-BON", "CRVAL1" => 20.0, "CRVAL2" => 120.0, "PV1_1" => 45.0))
    generic = parse_fits_wcs(generic_header)
    generic_world = pixel_to_world(generic, 512.0, 512.0)
    @test generic_world ≈ SVector(20.0, 120.0) atol=2e-12
    @test world_to_pixel(generic, generic_world) ≈ SVector(512.0, 512.0) atol=1e-10
    @test world_axis_types(generic) == (SPATIAL, SPATIAL)
    @test world_axis_physical_types(generic) == ("custom:xlat", "custom:xlon")
    @test output_frame(generic).reference_frame.radesys == ""
    @test_throws "not a matching longitude/latitude pair" parse_fits_wcs(merge(tan_header, Dict("CTYPE1" => "XLON-TAN", "CTYPE2" => "YLAT-TAN")))
    @test_throws "must be paired" parse_fits_wcs(merge(tan_header, Dict("CTYPE1" => "XXLN-TAN", "CTYPE2" => "Y")))

    galactic_header = merge(tan_header, Dict("CTYPE1" => "GLON-CAR", "CTYPE2" => "GLAT-CAR", "RADESYS" => "FK5", "EQUINOX" => 2000.0))
    galactic = output_frame(parse_fits_wcs(galactic_header)).reference_frame
    @test galactic.radesys == ""
    @test galactic.equinox === nothing
    linear_reference = output_frame(parse_fits_wcs(merge(linear_header, Dict("RADESYS" => "FK5", "EQUINOX" => 2000.0)))).reference_frame
    @test linear_reference.radesys == ""
    @test linear_reference.equinox === nothing

    log_reference = 10.0
    log_offset = 2.0
    generic_log_header = Dict("NAXIS" => 1, "CTYPE1" => "X----LOG", "CRPIX1" => 1.0, "CRVAL1" => log_reference, "CDELT1" => 2.0)
    generic_log = parse_fits_wcs(generic_log_header)
    generic_log_world = log_reference * exp(log_offset / log_reference)
    @test pixel_to_world(generic_log, 2.0) ≈ SVector(generic_log_world)
    @test world_to_pixel(generic_log, generic_log_world) ≈ SVector(2.0)
    @test all(isnan, world_to_pixel(generic_log, 0.0))
    @test_throws DomainError world_to_pixel(generic_log, 0.0; strict=true)
    @test_throws "reference value must be finite and positive" parse_fits_wcs(merge(generic_log_header, Dict("CRVAL1" => 0.0)))
    linear_negative = parse_fits_wcs(merge(generic_log_header, Dict("CTYPE1" => "X", "CRVAL1" => -10.0)))
    @test pixel_to_world(linear_negative, 1.0) == SVector(-10.0)

    time_header = Dict(
        "NAXIS" => 1, "CTYPE1" => "TIME", "TIMEUNIT" => "d", "TIMESYS" => "UTC", "TREFPOS" => "TOPOCENTER",
        "MJDREFI" => 59000, "MJDREFF" => 0.25, "OBSGEO-X" => 1.0, "OBSGEO-Y" => 2.0, "OBSGEO-Z" => 3.0,
        "CRPIX1" => 1.0, "CRVAL1" => 1.0, "CDELT1" => 0.5,
    )
    time = parse_fits_wcs(time_header)
    time_log = parse_fits_wcs(merge(time_header, Dict("CTYPE1" => "TIME-LOG")))
    @test pixel_to_world(time_log, 1.0) ≈ SVector(86400.0)
    @test pixel_to_world(time_log, 2.0) ≈ SVector(86400.0 * exp(0.5))
    @test world_to_pixel(time_log, 86400.0 * exp(0.5)) ≈ SVector(2.0)
    time_metadata = output_frame(time).reference_frame.time
    @test pixel_to_world(time, 1.0) == SVector(86400.0)
    @test pixel_to_world(parse_fits_wcs(time_header; preserve_units=true), 2.0) == SVector(1.5)
    @test world_axis_units(time) == ("s",)
    @test world_axis_units(parse_fits_wcs(time_header; preserve_units=true)) == ("d",)
    # Julian year and century use WCSLIB's factors; CUNITi takes precedence over TIMEUNIT.
    @test pixel_to_world(parse_fits_wcs(merge(time_header, Dict("TIMEUNIT" => "a"))), 2.0) == SVector(1.5 * 31557600.0)
    @test pixel_to_world(parse_fits_wcs(merge(time_header, Dict("TIMEUNIT" => "yr"))), 2.0) == SVector(1.5 * 31557600.0)
    @test pixel_to_world(parse_fits_wcs(merge(time_header, Dict("CUNIT1" => "cy"))), 2.0) == SVector(1.5 * 3155760000.0)
    @test world_axis_units(parse_fits_wcs(merge(time_header, Dict("CUNIT1" => "cy")); preserve_units=true)) == ("cy",)
    @test_throws "unsupported time CUNIT1" parse_fits_wcs(merge(time_header, Dict("CUNIT1" => "fortnight")))
    @test time_metadata.reference_epoch == 59000.25
    @test time_metadata.time_scale == "UTC"
    @test time_metadata.reference_position == "TOPOCENTER"
    @test time_metadata.location == (1.0, 2.0, 3.0)

    spectral_metadata_header = Dict("NAXIS" => 1, "CTYPE1" => "FREQ", "RESTFRQ" => 1.420405751e9)
    @test output_frame(parse_fits_wcs(spectral_metadata_header)).reference_frame.restfrq == 1.420405751e9

    spectral_cases = (
        ("FREQ-LOG", "Hz", 1.4e9, 1e6, (), 1.4015008038584952e9),
        ("WAVE-F2W", "m", 0.21, 1e-4, (), 0.21015010721944244),
        ("FREQ-W2F", "Hz", 1.4e9, 1e6, (), 1.4015016088666427e9),
        ("VELO-F2V", "m/s", 1000.0, 100.0, ("RESTFRQ" => 1.420405751e9,), 1150.000037502141),
        ("FREQ-V2F", "Hz", 1.4e9, 1e6, ("RESTFRQ" => 1.420405751e9,), 1.4015007811542246e9),
        ("VELO-W2V", "m/s", 1000.0, 100.0, ("RESTWAV" => 0.21,), 1149.999962408919),
        ("WAVE-V2W", "m", 0.21, 1e-4, ("RESTWAV" => 0.21,), 0.2101500536097144),
        ("AWAV-F2A", "m", 5e-7, 1e-10, (), 5.00150045012569e-7),
        ("FREQ-A2F", "Hz", 6e14, 1e10, (), 6.000150003750065e14),
        ("AWAV-W2A", "m", 5e-7, 1e-10, (), 5.001499999996591e-7),
        ("WAVE-A2W", "m", 5e-7, 1e-10, (), 5.001500000003412e-7),
        ("AWAV-V2A", "m", 5e-7, 1e-10, ("RESTFRQ" => 6e14,), 5.001500225498197e-7),
        ("VELO-A2V", "m/s", 1000.0, 100.0, ("RESTFRQ" => 6e14,), 1149.9999624797576),
        ("AFRQ-W2F", "rad/s", 2pi * 1.4e9, 2pi * 1e6, (), 8.80589431681944e9),
        ("ENER-W2F", "J", 6.62607015e-34 * 1.4e9, 6.62607015e-34 * 1e6, (), 9.286447975688238e-25),
        ("WAVN-W2F", "m-1", 1.4e9 / 299792458, 1e6 / 299792458, (), 4.674906160803561),
        ("VRAD-W2F", "m/s", 1000.0, 100.0, ("RESTFRQ" => 1.420405751e9,), 1149.9999249680382),
        ("VOPT-F2W", "m/s", 1000.0, 100.0, ("RESTWAV" => 0.21,), 1150.0000750422478),
        ("ZOPT-F2W", "", 0.01, 0.001, ("RESTWAV" => 0.21,), 0.0115022310361923),
        ("BETA-F2V", "", 0.001, 0.0001, ("RESTFRQ" => 1.420405751e9,), 0.0011500112275077232),
    )
    for (ctype_value, unit, crval_value, cdelt_value, extra, expected) in spectral_cases
        header = Dict{String,Any}(
            "NAXIS" => 1, "CTYPE1" => ctype_value, "CUNIT1" => unit,
            "CRPIX1" => 1.0, "CRVAL1" => crval_value, "CDELT1" => cdelt_value,
        )
        foreach(pair -> setindex!(header, pair.second, pair.first), extra)
        spectral = parse_fits_wcs(header)
        value = pixel_to_world(spectral, 2.5)
        @test only(value) ≈ expected rtol=8e-11
        @test world_to_pixel(spectral, value) ≈ SVector(2.5) atol=4e-9
    end
    log_ghz_header = Dict("NAXIS" => 1, "CTYPE1" => "FREQ-LOG", "CUNIT1" => "GHz", "CRPIX1" => 1.0, "CRVAL1" => 1.4, "CDELT1" => 0.001)
    log_ghz = parse_fits_wcs(log_ghz_header)
    log_ghz_preserved = parse_fits_wcs(log_ghz_header; preserve_units=true)
    @test pixel_to_world(log_ghz, 2.0) ≈ SVector(1.4010003572279062e9)
    @test pixel_to_world(log_ghz_preserved, 2.0) ≈ SVector(1.4010003572279062)
    @test world_axis_units(log_ghz) == ("Hz",)
    @test world_axis_units(log_ghz_preserved) == ("GHz",)
    @test @inferred(pixel_to_world(log_ghz, SVector(2.0))) ≈ SVector(1.4010003572279062e9)
    pixel_to_world(log_ghz, SVector(2.0))
    world_to_pixel(log_ghz, SVector(1.4010003572279062e9))
    @test allocated_fits_forward(log_ghz, SVector(2.0)) == 0
    @test allocated_fits_inverse(log_ghz, SVector(1.4010003572279062e9)) == 0
    @test isnan(only(pixel_to_world(log_ghz, NaN)))
    @test isnan(only(world_to_pixel(log_ghz, NaN)))
    @test_throws "incompatible with algorithm" parse_fits_wcs(Dict("NAXIS" => 1, "CTYPE1" => "FREQ-F2W", "CRVAL1" => 1.0))
    @test_throws "requires RESTFRQ or RESTWAV" parse_fits_wcs(Dict("NAXIS" => 1, "CTYPE1" => "VELO-F2V", "CRVAL1" => 0.0))

    grism_parameters = Dict{String,Any}("PV1_0" => 3.16e5, "PV1_1" => 1.0, "PV1_2" => 13.9)
    grism_cases = (
        (Dict{String,Any}("CTYPE1" => "WAVE-GRI", "CRVAL1" => 5.2252e-7, "CDELT1" => -4.334e-11, "CUNIT1" => "m"), 5.224766600224223e-7),
        (Dict{String,Any}("CTYPE1" => "FREQ-GRI", "CRVAL1" => 5.737e14, "CDELT1" => -4.758e7, "CUNIT1" => "Hz", "RESTFRQ" => 5.737e14), 5.736999524200039e14),
        (Dict{String,Any}("CTYPE1" => "VELO-GRI", "CRVAL1" => 0.0, "CDELT1" => 1000.0, "CUNIT1" => "m/s", "RESTFRQ" => 5.737e14), 999.9983529791468),
        (Dict{String,Any}("CTYPE1" => "VELO-GRA", "CRVAL1" => 0.0, "CDELT1" => 1000.0, "CUNIT1" => "m/s", "RESTFRQ" => 5.737e14), 999.9983530345745),
    )
    for (specific, expected) in grism_cases
        grism = parse_fits_wcs(merge(Dict{String,Any}("NAXIS" => 1, "CRPIX1" => 1.0), specific, grism_parameters))
        value = pixel_to_world(grism, 2.0)
        @test only(value) ≈ expected rtol=8e-11
        @test world_to_pixel(grism, value) ≈ SVector(2.0) atol=4e-9
    end
    grism_angstrom_header = Dict(
        "NAXIS" => 1, "CTYPE1" => "AWAV-GRA", "CRPIX1" => 1801.7, "CRVAL1" => 5225.2, "CDELT1" => -0.4334,
        "CUNIT1" => "Angstrom", "PV1_0" => 3.16e5, "PV1_1" => 1.0, "PV1_2" => 13.9,
    )
    grism_angstrom = parse_fits_wcs(grism_angstrom_header)
    @test pixel_to_world(grism_angstrom, 2000.0) ≈ SVector(5.139265914152931e-7) rtol=2e-15
    @test pixel_to_world(parse_fits_wcs(grism_angstrom_header; preserve_units=true), 2000.0) ≈ SVector(5139.265914152931) rtol=2e-15
    pixel_to_world(grism_angstrom, SVector(2000.0))
    world_to_pixel(grism_angstrom, SVector(5.139265914152931e-7))
    @test allocated_fits_forward(grism_angstrom, SVector(2000.0)) == 0
    @test allocated_fits_inverse(grism_angstrom, SVector(5.139265914152931e-7)) == 0
    @test isnan(only(pixel_to_world(grism_angstrom, NaN)))
    @test isnan(only(world_to_pixel(grism_angstrom, NaN)))
    @test_throws "singular wavelength equation" parse_fits_wcs(Dict("NAXIS" => 1, "CTYPE1" => "WAVE-GRI", "CRVAL1" => 5e-7))

    mixed_header = Dict(
        "NAXIS" => 4, "CTYPE1" => "FREQ-LOG", "CTYPE2" => "DEC--CAR", "CTYPE3" => "TIME", "CTYPE4" => "RA---CAR",
        "CUNIT1" => "GHz", "CUNIT2" => "deg", "CUNIT3" => "d", "CUNIT4" => "deg",
        "CRPIX1" => 1.0, "CRPIX2" => 1.0, "CRPIX3" => 1.0, "CRPIX4" => 1.0,
        "CRVAL1" => 1.4, "CRVAL2" => -20.0, "CRVAL3" => 2.0, "CRVAL4" => 30.0,
        "CDELT1" => 0.001, "CDELT2" => 0.1, "CDELT3" => 0.25, "CDELT4" => -0.1,
    )
    mixed = parse_fits_wcs(mixed_header)
    @test pixel_to_world(mixed, 2.0, 1.0, 3.0, 1.0) ≈ SVector(1.4010003572279062e9, -20.0, 216000.0, 30.0) atol=2e-7
    @test world_to_pixel(mixed, pixel_to_world(mixed, 2.0, 1.0, 3.0, 1.0)) ≈ SVector(2.0, 1.0, 3.0, 1.0) atol=2e-10

    @test pixel_to_world(parse_fits_wcs(Dict("" => nothing, "NAXIS" => 1, "CTYPE1" => "X", "CRVAL1" => 1.0)), 0.0) == SVector(1.0)
    fitsio_header = FITSIO.FITSHeader(["NAXIS", "", "CTYPE1", "CRVAL1"], Any[1, nothing, "X", 2.0], ["", "", "", ""])
    @test pixel_to_world(parse_fits_wcs(fitsio_header), 0.0) == SVector(2.0)
    fitsfiles_cards = FITSFiles.Card[FITSFiles.Card("NAXIS", 1), FITSFiles.Card(), FITSFiles.Card("CTYPE1", "X"), FITSFiles.Card("CRVAL1", 3.0)]
    @test pixel_to_world(parse_fits_wcs(fitsfiles_cards), 0.0) == SVector(3.0)

    tpd_card_keys = ["NAXIS", "CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2", "CPDIS1", "CPDIS2",
        "DP1", "DP1", "DP1", "DP1", "DP2", "DP2", "DP2", "DP2"]
    tpd_card_values = Any[2, "LINEAR", "LINEAR", 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, "TPD", "TPD",
        "NAXES: 2", "AXIS.1: 1", "AXIS.2: 2", "TPD.FWD.4: 0.01",
        "NAXES: 2", "AXIS.1: 2", "AXIS.2: 1", "TPD.FWD.4: -0.02"]
    tpd_fitsfiles_cards = FITSFiles.Card[FITSFiles.Card(key, value) for (key, value) in zip(tpd_card_keys, tpd_card_values)]
    @test pixel_to_world(parse_fits_wcs(tpd_fitsfiles_cards), 2.0, 3.0) ≈ SVector(2.04, 2.82) atol=2e-15
    tpd_fitsio_header = FITSIO.FITSHeader(tpd_card_keys, tpd_card_values, fill("", length(tpd_card_keys)))
    @test pixel_to_world(parse_fits_wcs(tpd_fitsio_header), 2.0, 3.0) ≈ SVector(2.04, 2.82) atol=2e-15

    sip_card_keys = ["NAXIS", "CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2", "CPDIS1", "CPDIS2",
        "DP1", "DP1", "DP1", "DP2", "DP2", "DP2"]
    sip_card_values = Any[2, "RA---TAN", "DEC--TAN", 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, "SIP", "SIP",
        "NAXES: 2", "SIP.FWD.2_0: 0.001", "SIP.REV.2_0: -0.001",
        "NAXES: 2", "SIP.FWD.0_2: 0.002", "SIP.REV.0_2: -0.002"]
    sip_expected = SVector(2.0031834025537707, 3.0133741951373705)
    sip_fitsfiles_cards = reduce(vcat, FITSFiles.Card(key, value) for (key, value) in zip(sip_card_keys, sip_card_values))
    @test pixel_to_world(parse_fits_wcs(sip_fitsfiles_cards), 2.0, 3.0) ≈ sip_expected atol=5e-14
    sip_fitsio_header = FITSIO.FITSHeader(sip_card_keys, sip_card_values, fill("", length(sip_card_keys)))
    @test pixel_to_world(parse_fits_wcs(sip_fitsio_header), 2.0, 3.0) ≈ sip_expected atol=5e-14

    lookup_cards = FITSFiles.Card[
        FITSFiles.Card("CTYPE1", "RA---TAN"), FITSFiles.Card("CTYPE2", "DEC--TAN"), FITSFiles.Card("CRPIX1", 2.0), FITSFiles.Card("CRPIX2", 2.0),
        FITSFiles.Card("CRVAL1", 30.0), FITSFiles.Card("CRVAL2", 10.0), FITSFiles.Card("CDELT1", -0.1), FITSFiles.Card("CDELT2", 0.1),
        FITSFiles.Card("D2IMDIS1", "LOOKUP"), FITSFiles.Card("D2IMERR1", 0.1), FITSFiles.Card("D2IM1", "NAXES: 2"),
        FITSFiles.Card("D2IM1", "AXIS.1: 1"), FITSFiles.Card("D2IM1", "AXIS.2: 2"),
    ]
    lookup_primary = FITSFiles.HDU(zeros(4, 4), lookup_cards)
    lookup_extension = FITSFiles.HDU(ones(1, 1), FITSFiles.Card[FITSFiles.Card("EXTNAME", "D2IMARR"), FITSFiles.Card("EXTVER", 1), FITSFiles.Card("CRPIX1", 0.0), FITSFiles.Card("CRPIX2", 0.0)])
    fitsfiles_lookup = parse_fits_wcs([lookup_primary, lookup_extension], 1)
    for distortion in ("TPD", "SIP")
        unsupported_cards = FITSFiles.Card[card.key == "D2IMDIS1" ? FITSFiles.Card("D2IMDIS1", distortion) : card for card in lookup_cards]
        @test_throws "unsupported Paper IV distortion type \"$distortion\" in D2IMDIS1" parse_fits_wcs([FITSFiles.HDU(zeros(4, 4), unsupported_cards), lookup_extension], 1)
    end
    fitsfiles_plain = parse_fits_wcs(merge(tan_header, Dict("CRPIX1" => 2.0, "CRPIX2" => 2.0, "CRVAL1" => 30.0, "CRVAL2" => 10.0, "CDELT1" => -0.1, "CDELT2" => 0.1)))
    @test pixel_to_world(fitsfiles_lookup, 2.0, 2.0) ≈ pixel_to_world(fitsfiles_plain, 3.0, 2.0) atol=1e-12
    lookup_1d_cards = FITSFiles.Card[card for card in lookup_cards if card.key != "D2IM1" || card.value != "NAXES: 2"]
    filter!(card -> card.key != "D2IM1" || card.value != "AXIS.2: 2", lookup_1d_cards)
    push!(lookup_1d_cards, FITSFiles.Card("D2IM1", "NAXES: 1"))
    lookup_1d_primary = FITSFiles.HDU(zeros(4, 4), lookup_1d_cards)
    lookup_1d_extension = FITSFiles.HDU([0.0, 1.0, 2.0], FITSFiles.Card[
        FITSFiles.Card("EXTNAME", "D2IMARR"), FITSFiles.Card("EXTVER", 1), FITSFiles.Card("CRPIX1", 1.0), FITSFiles.Card("CRVAL1", 1.0), FITSFiles.Card("CDELT1", 1.0),
    ])
    fitsfiles_lookup_1d = parse_fits_wcs([lookup_1d_primary, lookup_1d_extension], 1)
    @test pixel_to_world(fitsfiles_lookup_1d, 2.0, 2.0) ≈ pixel_to_world(fitsfiles_plain, 3.0, 2.0) atol=1e-12
    @test pixel_to_world(parse_fits_wcs([lookup_primary, lookup_extension], 1; minerr=0.2), 2.0, 2.0) ≈ pixel_to_world(fitsfiles_plain, 2.0, 2.0) atol=1e-12
    malformed_cards = FITSFiles.Card[FITSFiles.Card("CTYPE1", "X"), FITSFiles.Card("WCSAXESA", 1), FITSFiles.Card("CTYPE1A", "DEC--TAN")]
    @test_throws "FITS WCS key 'A' is malformed" find_all_fits_wcs([FITSFiles.HDU(zeros(2), malformed_cards)], 1)
    wcsaxes_only = Dict("WCSAXES" => 1, "WCSAXESB" => 2)
    @test keys(find_all_fits_wcs(wcsaxes_only)) == Set([' ', 'B'])
    @test input_dimension(find_all_fits_wcs(wcsaxes_only)[' '].forward) == 1
    @test input_dimension(find_all_fits_wcs(wcsaxes_only)['B'].forward) == 2

    tab_primary = FITSFiles.HDU(zeros(3), FITSFiles.Card[
        FITSFiles.Card("CTYPE1", "FREQ-TAB"), FITSFiles.Card("CUNIT1", "GHz"), FITSFiles.Card("CRPIX1", 1.0),
        FITSFiles.Card("CRVAL1", 1.0), FITSFiles.Card("CDELT1", 1.0), FITSFiles.Card("PS1_0", "WCS-TABLE"), FITSFiles.Card("PS1_1", "coord"),
    ])
    tab_extension = FITSFiles.HDU((COORD=[10.0, 20.0, 40.0],), FITSFiles.Card[
        FITSFiles.Card("EXTNAME", "WCS-TABLE"), FITSFiles.Card("EXTVER", 1), FITSFiles.Card("EXTLEVEL", 1),
    ])
    fitsfiles_tab = parse_fits_wcs([tab_primary, tab_extension], 1)
    padded_tab_primary = FITSFiles.HDU(tab_primary.data, FITSFiles.Card[card.key == "PS1_0" ? FITSFiles.Card("PS1_0", "WCS-TABLE ") : card for card in tab_primary.cards])
    @test pixel_to_world(parse_fits_wcs([padded_tab_primary, tab_extension], 1), 2.5) == SVector(3.0e10)
    @test pixel_to_world(fitsfiles_tab, 2.5) == SVector(3.0e10)
    @test world_to_pixel(fitsfiles_tab, 3.0e10) == SVector(2.5)

    # Repeated index values define no interval, and mixed-dimensional tables must stay allocation-free.
    repeated_index = FITSTabularTable((1,), reshape([1.0, 4.0, 9.0, 16.0], 1, :), ([1.0, 2.0, 2.0, 3.0],))
    decreasing_index = FITSTabularTable((1,), reshape([1.0, 4.0, 9.0, 16.0], 1, :), ([3.0, 2.0, 2.0, 1.0],))
    @test FITSTabularTransform(Val(1), (repeated_index,), (1.0,))(SVector(2.0)) == SVector(4.0)
    @test FITSTabularTransform(Val(1), (repeated_index,), (1.0,))(SVector(2.5)) == SVector(12.5)
    @test FITSTabularTransform(Val(1), (decreasing_index,), (1.0,))(SVector(2.0)) == SVector(4.0)
    coupled_coordinates = [axis == 1 ? i + 0.1j : j + 0.1i for axis in 1:2, i in 1:5, j in 1:6]
    coupled_table = FITSTabularTable((2, 3), coupled_coordinates, (collect(1.0:5.0), collect(1.0:6.0)))
    mixed_tables = FITSTabularTransform(Val(3), (repeated_index, coupled_table), (1.0, 1.0, 1.0))
    @test mixed_tables(SVector(2.5, 2.0, 3.0)) ≈ SVector(12.5, 2.3, 3.2)
    allocated_tables(transform, value) = @allocated transform(value)
    allocated_tables(mixed_tables, SVector(2.5, 2.0, 3.0))
    @test allocated_tables(mixed_tables, SVector(2.5, 2.0, 3.0)) == 0

    mktempdir(@__DIR__) do directory
        read_container_wcs(path; key=' ') = (
            FITSIO.FITS(path) do file
                parse_fits_wcs(file[1]; key)
            end,
            parse_fits_wcs(FITSFiles.fits(path), 1; key),
        )

        path = joinpath(directory, "lookup.fits")
        lookup_keywords = ["NAXIS", "CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2",
            "D2IMDIS1", "D2IMDIS2", "D2IM1", "D2IM1", "D2IM1", "D2IM1", "D2IM2", "D2IM2", "D2IM2", "D2IM2", "CPDIS1", "CPDIS2", "DP1", "DP1", "DP1", "DP1", "DP2", "DP2", "DP2", "DP2",
            "WCSAXES", "CTYPE1A", "CTYPE2A", "CRPIX1A", "CRPIX2A", "CRVAL1A", "CRVAL2A", "CDELT1A", "CDELT2A", "D2IMDIS2A", "D2IM2A", "D2IM2A", "D2IM2A", "D2IM2A"]
        lookup_values = Any[2, "LINEAR", "LINEAR", 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
            "LOOKUP", "LOOKUP", "NAXES: 2", "AXIS.1: 1", "AXIS.2: 2", "EXTVER: 1", "NAXES: 2", "AXIS.1: 2", "AXIS.2: 1", "EXTVER: 2", "LOOKUP", "LOOKUP",
            "NAXES: 2", "AXIS.1: 1", "AXIS.2: 2", "EXTVER: 3", "NAXES: 2", "AXIS.1: 2", "AXIS.2: 1", "EXTVER: 4",
            2, "LINEAR", "LINEAR", 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, "LOOKUP", "NAXES: 2", "AXIS.1: 2", "AXIS.2: 1", "EXTVER: 9"]
        lookup_header = FITSIO.FITSHeader(lookup_keywords[1:9], lookup_values[1:9], fill("", 9))
        FITSIO.FITS(path, "w") do file
            FITSIO.write(file, zeros(2, 2); header=lookup_header)
            for (keyword, value) in zip(lookup_keywords[10:end], lookup_values[10:end])
                card_value = value isa AbstractString ? "'$(value)'" : string(value)
                FITSIO.Libcfitsio.fits_write_record(getfield(file[1], :fitsfile), rpad(keyword, 8) * "= $card_value")
            end
            for (name, version, value) in (("D2IMARR", 1, 1.0), ("D2IMARR", 2, 10.0), ("WCSDVARR", 3, 100.0), ("WCSDVARR", 4, 1000.0), ("D2IMARR", 8, -999.0))
                FITSIO.write(file, fill(value, 2, 2); name, ver=version)
            end
            FITSIO.write(file, [10.0 20.0; 30.0 40.0]; name="D2IMARR", ver=9)
        end
        FITSIO.FITS(path) do file
            primary_wcs = parse_fits_wcs(file[1])
            alternate_wcs = parse_fits_wcs(file[1]; key='A')
            @test pixel_to_world(primary_wcs, 2.0, 3.0) == SVector(103.0, 1013.0)
            @test pixel_to_world(alternate_wcs, 1.0, 2.0) == SVector(1.0, 32.0)
        end

        lookup_1d_path = joinpath(directory, "lookup_1d.fits")
        lookup_1d_header = FITSIO.FITSHeader(["NAXIS", "CTYPE1", "CTYPE2", "D2IMDIS1"], Any[2, "LINEAR", "LINEAR", "LOOKUP"], fill("", 4))
        FITSIO.FITS(lookup_1d_path, "w") do file
            FITSIO.write(file, zeros(2, 2); header=lookup_1d_header)
            for value in ("NAXES: 1", "AXIS.1: 1")
                FITSIO.Libcfitsio.fits_write_record(getfield(file[1], :fitsfile), "D2IM1   = '$(value)'")
            end
            FITSIO.write(file, [0.0, 1.0, 2.0]; name="D2IMARR", ver=1)
        end
        for wcs in read_container_wcs(lookup_1d_path)
            @test pixel_to_world(wcs, 2.0, 4.0) == SVector(3.0, 4.0)
        end

        duplicate_path = joinpath(directory, "duplicate.fits")
        duplicate_header = FITSIO.FITSHeader(["NAXIS", "CTYPE1", "CTYPE2", "D2IMDIS1"], Any[2, "LINEAR", "LINEAR", "LOOKUP"], fill("", 4))
        FITSIO.FITS(duplicate_path, "w") do file
            FITSIO.write(file, zeros(2, 2); header=duplicate_header)
            for value in ("NAXES: 2", "AXIS.1: 1", "AXIS.2: 2")
                FITSIO.Libcfitsio.fits_write_record(getfield(file[1], :fitsfile), "D2IM1   = '$(value)'")
            end
            FITSIO.write(file, ones(2, 2); name="D2IMARR", ver=1)
            FITSIO.write(file, zeros(2, 2); name="D2IMARR", ver=1)
        end
        FITSIO.FITS(duplicate_path) do file
            @test_throws "is ambiguous (2 matching HDUs)" parse_fits_wcs(file[1])
        end
        @test_throws "is ambiguous (2 matching HDUs)" parse_fits_wcs(FITSFiles.fits(duplicate_path), 1)

        tab_path = joinpath(directory, "tab.fits")
        tab_image_header = FITSIO.FITSHeader(
            ["NAXIS", "CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2", "PS1_0", "PS2_0", "PS1_1", "PS2_1", "PV1_3", "PV2_3"],
            Any[2, "RA---TAB", "DEC--TAB", 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, "WCS-TABLE", "WCS-TABLE", "COORD", "COORD", 1, 2], fill("", 15),
        )
        coordinates = reshape([11.0, 11.0, 12.0, 21.0, 21.0, 12.0, 22.0, 22.0], 8, 1)
        tab_table_header = FITSIO.FITSHeader(["TDIM1", "EXTLEVEL"], Any["(2,2,2)", 1], ["", ""])
        FITSIO.FITS(tab_path, "w") do file
            FITSIO.write(file, zeros(2, 2); header=tab_image_header)
            FITSIO.write(file, Dict("COORD" => coordinates); header=tab_table_header, name="WCS-TABLE", ver=1)
        end
        for wcs in read_container_wcs(tab_path)
            @test pixel_to_world(wcs, 1.5, 1.5) == SVector(16.5, 16.5)
            @test world_to_pixel(wcs, 16.5, 16.5; initial_guess=(1.5, 1.5)) ≈ SVector(1.5, 1.5) atol=1e-10
        end

        # FITSFiles parses short real literals such as 150.1 as Float32; the adapter must recover the decimal value.
        short_real_path = joinpath(directory, "short_reals.fits")
        FITSIO.FITS(short_real_path, "w") do file
            FITSIO.write(file, zeros(2, 2); header=FITSIO.FITSHeader(["CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CD1_1", "CD2_2"],
                Any["RA---TAN", "DEC--TAN", 50.5, 40.5, 150.1, 2.2, -1.1e-4, 1.05e-4], fill("", 8)))
        end
        short_real_header = Dict("NAXIS" => 2, "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CRPIX1" => 50.5, "CRPIX2" => 40.5, "CRVAL1" => 150.1, "CRVAL2" => 2.2, "CD1_1" => -1.1e-4, "CD2_2" => 1.05e-4)
        for wcs in read_container_wcs(short_real_path)
            @test pixel_to_world(wcs, 10.0, 20.0) ≈ pixel_to_world(parse_fits_wcs(short_real_header), 10.0, 20.0) atol=1e-12
        end

        # Header-only CPDIS SIP records must not be treated as linked lookup tables.
        sip_path = joinpath(directory, "sip_records.fits")
        FITSIO.FITS(sip_path, "w") do file
            FITSIO.write(file, zeros(2, 2); header=FITSIO.FITSHeader(sip_card_keys[2:11], sip_card_values[2:11], fill("", 10)))
            for (keyword, value) in zip(sip_card_keys[12:end], sip_card_values[12:end])
                FITSIO.Libcfitsio.fits_write_record(getfield(file[1], :fitsfile), rpad(keyword, 8) * "= '$(value)'")
            end
        end
        for wcs in read_container_wcs(sip_path)
            @test pixel_to_world(wcs, 2.0, 3.0) ≈ sip_expected atol=5e-14
        end
        @test pixel_to_world(find_all_fits_wcs(FITSFiles.fits(sip_path), 1)[' '], 2.0, 3.0) ≈ sip_expected atol=5e-14
        FITSIO.FITS(sip_path) do file
            @test pixel_to_world(find_all_fits_wcs(file[1])[' '], 2.0, 3.0) ≈ sip_expected atol=5e-14
        end
    end

    alternate_header = merge(linear_header, Dict(
        "WCSAXESA" => 2, "CTYPE1A" => "X", "CTYPE2A" => "Y",
        "CRPIX1A" => 1.0, "CRPIX2A" => 1.0, "CRVAL1A" => 10.0, "CRVAL2A" => 20.0,
        "CDELT1A" => 1.0, "CDELT2A" => 1.0, "WCSNAMEA" => "alternate",
    ))
    alternate = parse_fits_wcs(alternate_header; key='A')
    @test alternate.name == "alternate"
    @test pixel_to_world(alternate, 1.0, 1.0) == SVector(10.0, 20.0)
    all_wcs = find_all_fits_wcs(alternate_header)
    @test Set(keys(all_wcs)) == Set((' ', 'A'))
    @test pixel_to_world(all_wcs[' '], 10.0, 20.0) == SVector(100.0, 200.0)
    @test_throws "FITS WCS key 'A' is malformed" find_all_fits_wcs(merge(alternate_header, Dict("CTYPE2A" => "DEC--TAN")))

    partial = parse_fits_wcs(Dict("NAXIS" => 2, "WCSAXES" => 3, "NAXIS1" => 10, "NAXIS2" => 20))
    @test partial.pixel_shape == (10, 20, nothing)
    @test WCS(Step(CoordinateFrame{2}(name="input"), IdentityTransform{2}()), Step(CoordinateFrame{2}(name="output")); pixel_shape=(10, nothing)).pixel_shape == (10, nothing)

    @test_throws "A_ORDER and B_ORDER are absent" parse_fits_wcs(merge(tan_header, Dict("CTYPE1" => "RA---TAN-SIP", "CTYPE2" => "DEC--TAN-SIP")))
    @test_throws "external distortion or table data" parse_fits_wcs(merge(linear_header, Dict("CPDIS1" => "LOOKUP")))
    @test_throws "unsupported nonlinear algorithm" parse_fits_wcs(Dict("NAXIS" => 1, "CTYPE1" => "FREQ-XYZ"))
    @test_throws "key must be" parse_fits_wcs(linear_header; key='a')

    sip_header = merge(tan_header, Dict(
        "CTYPE1" => "RA---TAN-SIP", "CTYPE2" => "DEC--TAN-SIP",
        "A_ORDER" => 2, "A_2_0" => 5e-6, "A_0_2" => 2e-6,
        "B_ORDER" => 2, "B_2_0" => 1e-6, "B_1_1" => 3e-6,
    ))
    sip = parse_fits_wcs(sip_header)
    # AP/BP map CRVAL exactly to CRPIX, the TAN native pole, where the chain-rule Jacobian is NaN.
    exact_reference_sip = Dict(
        "NAXIS" => 2, "CTYPE1" => "RA---TAN-SIP", "CTYPE2" => "DEC--TAN-SIP", "CRPIX1" => 512.0, "CRPIX2" => 512.0,
        "CRVAL1" => 10.0, "CRVAL2" => 20.0, "CDELT1" => 1.0, "CDELT2" => 1.0, "A_ORDER" => 2, "B_ORDER" => 2, "A_2_0" => 0.001, "B_0_2" => 0.002,
        "AP_ORDER" => 2, "BP_ORDER" => 2, "AP_2_0" => -0.001, "BP_0_2" => -0.002,
    )
    @test world_to_pixel(parse_fits_wcs(exact_reference_sip), 10.0, 20.0) ≈ SVector(512.0, 512.0) atol=1e-10
    sip_world = pixel_to_world(sip, 400.0, 300.0)
    @test sip_world ≈ SVector(83.8533100335119, -5.449965264308421) atol=2e-11
    @test world_to_pixel(sip, sip_world; initial_guess=SVector(400.0, 300.0)) ≈ SVector(400.0, 300.0) atol=1e-8

    sip_inverse_header = merge(sip_header, Dict(
        "AP_ORDER" => 2, "AP_2_0" => -5e-6, "AP_0_2" => -2e-6,
        "BP_ORDER" => 2, "BP_2_0" => -1e-6, "BP_1_1" => -3e-6,
    ))
    sip_inverse = parse_fits_wcs(sip_inverse_header)
    @test !has_inverse(sip_inverse.forward)
    sip_inverse_world = pixel_to_world(sip_inverse, 400.0, 300.0)
    sip_reverse_only = GWCS._inverse_initialize(sip_inverse.forward, sip_inverse_world)[1]
    @test maximum(abs, pixel_to_world(sip_inverse, sip_reverse_only) - sip_inverse_world) > 1e-8
    @test world_to_pixel(sip_inverse, sip_inverse_world) ≈ SVector(400.0, 300.0) atol=1e-8

    # Strong corner distortion cross-checked directly with WCSLIB 8.9 wcsp2s/wcss2p.
    strong_sip_header = Dict(
        "NAXIS" => 2, "NAXIS1" => 512, "NAXIS2" => 512,
        "CTYPE1" => "RA---TAN-SIP", "CTYPE2" => "DEC--TAN-SIP",
        "CRPIX1" => 256.0, "CRPIX2" => 256.0, "CRVAL1" => 120.0, "CRVAL2" => 35.0,
        "CD1_1" => -0.01, "CD2_2" => 0.01,
        "A_ORDER" => 2, "A_2_0" => 5e-4, "B_ORDER" => 2, "B_0_2" => -5e-4,
        "AP_ORDER" => 3, "AP_2_0" => -5e-4, "AP_3_0" => 5e-7,
        "BP_ORDER" => 3, "BP_0_2" => 5e-4, "BP_0_3" => 5e-7,
    )
    strong_sip = parse_fits_wcs(strong_sip_header)
    strong_sip_pixel = SVector(1.0, 1.0)
    strong_sip_world = pixel_to_world(strong_sip, strong_sip_pixel)
    @test strong_sip_world ≈ SVector(122.62204455315441, 32.10026082474328) atol=3e-14
    strong_sip_initial = GWCS._inverse_initialize(strong_sip.forward, strong_sip_world)[1]
    @test maximum(abs, strong_sip_initial - strong_sip_pixel) > 3
    @test world_to_pixel(strong_sip, strong_sip_world) ≈ strong_sip_pixel atol=2e-10

    sip_cube_header = merge(sip_header, Dict("NAXIS" => 3, "CTYPE3" => "FREQ", "CUNIT3" => "Hz", "CRPIX3" => 1.0, "CRVAL3" => 10.0, "CDELT3" => 2.0))
    @test pixel_to_world(parse_fits_wcs(sip_cube_header), 512.0, 512.0, 2.0) ≈ SVector(83.8221, -5.3911, 12.0) atol=1e-11

    lookup_table = BilinearLookupTable(ones(1, 1); crpix=(0, 0))
    lookup_header = merge(tan_header, Dict("D2IMDIS1" => "LOOKUP", "D2IM1.NAXES" => 2, "D2IM1.AXIS.1" => 1, "D2IM1.AXIS.2" => 2))
    lookup_data = (detector_to_image=(lookup_table, nothing), cpdis=(nothing, nothing))
    lookup_wcs = GWCS._parse_fits_wcs(GWCS._fits_header(lookup_header); lookup_data)
    @test pixel_to_world(lookup_wcs, 512.0, 512.0) ≈ pixel_to_world(tan, 513.0, 512.0) atol=1e-12
    one_dimensional_tpd_header = Dict(
        "NAXIS" => 1, "CTYPE1" => "LINEAR", "CRPIX1" => 1.0, "CRVAL1" => 0.0, "CDELT1" => 1.0,
        "CPDIS1" => "TPD", "DP1.NAXES" => 1, "DP1.TPD.FWD.1" => 1.0, "D2IMDIS1" => "LOOKUP",
    )
    @test_throws "FITS lookup distortions require at least two pixel axes" GWCS._parse_fits_wcs(GWCS._fits_header(one_dimensional_tpd_header); lookup_data)
    pixel_to_world(lookup_wcs, SVector(512.0, 512.0))
    @test allocated_fits_forward(lookup_wcs, SVector(512.0, 512.0)) == 0
    @test all(isnan, pixel_to_world(lookup_wcs, NaN, 512.0))
    bilinear = BilinearLookupTable([0.0 1.0; 2.0 3.0])
    @test bilinear(SVector(0.5, 0.5)) == SVector(1.5)
    @test_throws "CDELT values must be nonzero" BilinearLookupTable(ones(2, 2); cdelt=(1, 0))

    pv_pole_header = Dict(
        "NAXIS" => 2, "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CRPIX1" => 1.0, "CRPIX2" => 1.0,
        "CRVAL1" => 10.0, "CRVAL2" => 20.0, "CDELT1" => 1.0, "CDELT2" => 1.0,
        "PV1_1" => 5.0, "PV1_2" => 30.0, "PV1_3" => 80.0, "PV1_4" => 70.0,
        "LONPOLE" => 10.0, "LATPOLE" => 60.0,
    )
    # WCSLIB gives PV1_3/PV1_4 precedence over LONPOLE/LATPOLE.
    pv_pole = parse_fits_wcs(pv_pole_header)
    @test pixel_to_world(pv_pole, 1.0, 1.0) ≈ SVector(307.10135177, 14.47667920) atol=6e-9
    @test pixel_to_world(pv_pole, 2.0, 3.0) ≈ SVector(304.88253948, 15.10335044) atol=7e-9
    pv_offset = parse_fits_wcs(merge(pv_pole_header, Dict("PV1_0" => 1.0)))
    @test pixel_to_world(pv_offset, 1.0, 1.0) ≈ SVector(10.0, 20.0) atol=2e-12
    @test pixel_to_world(pv_offset, 2.0, 3.0) ≈ SVector(9.54764649, 20.63728413) atol=6e-9

    ordered_header = merge(tan_header, Dict(
        "CRPIX1" => 1.0, "CRPIX2" => 1.0, "CTYPE1" => "RA---TAN-SIP", "CTYPE2" => "DEC--TAN-SIP",
        "A_ORDER" => 2, "A_2_0" => 0.1, "B_ORDER" => 2,
    ))
    ordered_lookup = (
        detector_to_image=(BilinearLookupTable(ones(1, 1)), nothing),
        cpdis=(BilinearLookupTable(fill(2.0, 1, 1)), nothing),
    )
    ordered = GWCS._parse_fits_wcs(GWCS._fits_header(ordered_header); lookup_data=ordered_lookup)
    @test get_transform(ordered, "detector", "intermediate")(SVector(2.0, 2.0)) ≈ SVector(-0.001222232, 0.00027778) atol=1e-15

    tab_header = Dict(
        "NAXIS" => 1, "CTYPE1" => "FREQ-TAB", "CUNIT1" => "GHz", "CRPIX1" => 1.0, "CRVAL1" => 1.0, "CDELT1" => 1.0,
        "PS1_0" => "WCS-TABLE", "PS1_1" => "COORD",
    )
    tabular = FITSTabularTable((1,), reshape([10.0, 20.0, 40.0], 1, 3), ([1.0, 2.0, 3.0],))
    tab_wcs = GWCS._parse_fits_wcs(GWCS._fits_header(tab_header); tabular_data=(tabular,))
    @test pixel_to_world(tab_wcs, 2.5) == SVector(3.0e10)
    @test world_to_pixel(tab_wcs, 3.0e10) == SVector(2.5)
    @test pixel_to_world(tab_wcs, 0.5) == SVector(5.0e9)
    @test world_to_pixel(tab_wcs, 5.0e9) == SVector(0.5)
    @test all(isnan, pixel_to_world(tab_wcs, 0.49))
    @test_throws DomainError pixel_to_world(tab_wcs, 0.49; strict=true)
    @test pixel_to_world(GWCS._parse_fits_wcs(GWCS._fits_header(tab_header); preserve_units=true, tabular_data=(tabular,)), 2.5) == SVector(30.0)
    pixel_to_world(tab_wcs, SVector(2.5))
    world_to_pixel(tab_wcs, SVector(3.0e10))
    @test allocated_fits_forward(tab_wcs, SVector(2.5)) == 0
    @test allocated_fits_inverse(tab_wcs, SVector(3.0e10)) == 0
    @test isnan(only(pixel_to_world(tab_wcs, NaN)))
    @test isnan(only(world_to_pixel(tab_wcs, NaN)))
    @test_throws "must be monotonic" FITSTabularTable((1,), reshape([10.0, 20.0, 40.0], 1, 3), ([1.0, 3.0, 2.0],))
    @test_throws "must be finite" FITSTabularTable((1,), reshape([10.0, NaN], 1, 2), ([1.0, 2.0],))

    coupled_header = Dict(
        "NAXIS" => 2, "CTYPE1" => "RA---TAB", "CTYPE2" => "DEC--TAB", "CRPIX1" => 1.0, "CRPIX2" => 1.0,
        "CRVAL1" => 1.0, "CRVAL2" => 1.0, "CDELT1" => 1.0, "CDELT2" => 1.0,
        "PS1_0" => "WCS-TABLE", "PS2_0" => "WCS-TABLE", "PS1_1" => "COORD", "PS2_1" => "COORD", "PV1_3" => 1, "PV2_3" => 2,
    )
    coupled_coordinates = zeros(2, 2, 2)
    for first in 1:2, second in 1:2
        coupled_coordinates[1, first, second] = 10first + second
        coupled_coordinates[2, first, second] = first + 10second
    end
    coupled_table = FITSTabularTable((1, 2), coupled_coordinates, ([1.0, 2.0], [1.0, 2.0]))
    coupled_wcs = GWCS._parse_fits_wcs(GWCS._fits_header(coupled_header); tabular_data=(coupled_table,))
    coupled_world = pixel_to_world(coupled_wcs, 1.5, 1.5)
    @test coupled_world == SVector(16.5, 16.5)
    @test world_to_pixel(coupled_wcs, coupled_world; initial_guess=(1.5, 1.5)) ≈ SVector(1.5, 1.5) atol=1e-10

    tpv_header = merge(tan_header, Dict(
        "CTYPE1" => "RA---TPV", "CTYPE2" => "DEC--TPV",
        "PV1_1" => 1.0, "PV1_4" => 0.01, "PV2_1" => 1.0, "PV2_4" => -0.02,
    ))
    tpv = parse_fits_wcs(tpv_header)
    tpv_world = pixel_to_world(tpv, 513.0, 514.0)
    @test tpv_world ≈ SVector(83.82181819590602, -5.390566662335204) atol=2e-11
    @test world_to_pixel(tpv, tpv_world; initial_guess=(513.0, 514.0)) ≈ SVector(513.0, 514.0) atol=1e-8
    # This TPV folds at pixel x = 462; like WCSLIB, the default start is the undistorted solution, which selects
    # the branch containing the original pixel rather than the mirror image reached from pixel (1, 1).
    @test world_to_pixel(tpv, pixel_to_world(tpv, 632.0, 407.5)) ≈ SVector(632.0, 407.5) atol=1e-8

    empty_tpv = parse_fits_wcs(filter(pair -> !startswith(pair.first, "PV"), tpv_header))
    @test pixel_to_world(empty_tpv, 513.0, 514.0) ≈ SVector(83.8221, -5.3911) atol=1e-12
    @test_throws "parameter numbers 0 through 39" parse_fits_wcs(merge(tpv_header, Dict("PV1_40" => 1.0)))
    @test_throws "Paper IV distortion function" parse_fits_wcs(merge(tan_header, Dict("CTYPE1" => "RA---TPD", "CTYPE2" => "DEC--TPD")))

    # WCSLIB applies sequent TPV distortion after PC but before CDELT.
    tpv_pc_header = Dict(
        "NAXIS" => 2, "CTYPE1" => "RA---TPV", "CTYPE2" => "DEC--TPV",
        "CRPIX1" => 128.0, "CRPIX2" => 96.0, "CRVAL1" => 120.0, "CRVAL2" => 35.0,
        "CDELT1" => -0.05, "CDELT2" => 0.05,
        "PC1_1" => 1.0, "PC1_2" => 0.1, "PC2_1" => -0.2, "PC2_2" => 1.0,
        "PV1_1" => 1.0, "PV1_4" => 0.01, "PV2_1" => 1.0, "PV2_4" => -0.02,
    )
    @test pixel_to_world(parse_fits_wcs(tpv_pc_header), 130.0, 99.0) ≈ SVector(119.85616564526542, 35.12315485197094) atol=3e-14

    # TPV polynomials use the longitude and latitude intermediate coordinates,
    # regardless of where those celestial axes occur in the FITS axis order.
    swapped_tpv_header = Dict(
        "NAXIS" => 2, "CTYPE1" => "DEC--TPV", "CTYPE2" => "RA---TPV",
        "CRPIX1" => 100.0, "CRPIX2" => 200.0, "CRVAL1" => -22.5, "CRVAL2" => 123.4,
        "CDELT1" => 0.035, "CDELT2" => -0.02,
        "PC1_1" => 1.1, "PC1_2" => -0.2, "PC2_1" => 0.3, "PC2_2" => 0.8,
        "PV1_0" => -0.02, "PV1_1" => 0.5, "PV1_2" => 0.5,
        "PV1_4" => 0.0004, "PV1_5" => -0.0003, "PV1_6" => 0.0004,
        "PV2_0" => 0.01, "PV2_1" => 1.02, "PV2_2" => -0.03,
        "PV2_4" => 0.0007, "PV2_5" => -0.0002, "PV2_6" => 0.0005,
    )
    swapped_models = GWCS._fits_tpv_models(GWCS._fits_header(swapped_tpv_header), 2, "", 2, 1)
    @test swapped_models[1].axes == (1, 2)
    @test swapped_models[2].axes == (2, 1)
    swapped_tpv = parse_fits_wcs(swapped_tpv_header)
    swapped_world = pixel_to_world(swapped_tpv, 112.3, 176.4)
    @test swapped_world ≈ SVector(-22.43599315380056, 123.73858159819254) atol=3e-14
    @test world_to_pixel(swapped_tpv, swapped_world; initial_guess=(112.3, 176.4)) ≈ SVector(112.3, 176.4) atol=1e-10

    cube_tpv_header = Dict(
        "NAXIS" => 3, "CTYPE1" => "FREQ", "CTYPE2" => "RA---TPV", "CTYPE3" => "DEC--TPV",
        "CRPIX1" => 10.0, "CRPIX2" => 100.0, "CRPIX3" => 200.0,
        "CRVAL1" => 1e9, "CRVAL2" => 123.4, "CRVAL3" => -22.5,
        "CDELT1" => 2e6, "CDELT2" => -0.02, "CDELT3" => 0.035,
        "PC1_1" => 1.0, "PC1_2" => 0.1, "PC1_3" => 0.0,
        "PC2_1" => 0.0, "PC2_2" => 0.8, "PC2_3" => 0.3,
        "PC3_1" => 0.0, "PC3_2" => -0.2, "PC3_3" => 1.1,
        "PV2_1" => 1.0, "PV2_4" => 0.0007, "PV3_1" => 1.0, "PV3_4" => -0.0004,
    )
    cube_models = GWCS._fits_tpv_models(GWCS._fits_header(cube_tpv_header), 3, "", 2, 3)
    @test cube_models[1] === nothing
    @test cube_models[2].axes == (2, 3)
    @test cube_models[3].axes == (3, 2)
    @test pixel_to_world(parse_fits_wcs(cube_tpv_header), 12.3, 112.3, 176.4) ≈ SVector(1.00706e9, 123.33969797345996, -23.50589278224176) atol=3e-14
    @test_throws "TPV cannot be combined" parse_fits_wcs(merge(tpv_header, Dict("CQDIS1" => "TPD", "DQ1.NAXES" => 1)))
    # TPV only couples the celestial axes; PC1_2 couples the spectral axis to pixel axis 2.
    @test axis_correlation_matrix(parse_fits_wcs(cube_tpv_header)) == SMatrix{3,3,Bool}([1 1 0; 0 1 1; 0 1 1])

    prior_tpd_header = Dict(
        "NAXIS" => 2, "CTYPE1" => "LINEAR", "CTYPE2" => "LINEAR",
        "CRPIX1" => 10.0, "CRPIX2" => 20.0, "CRVAL1" => 5.0, "CRVAL2" => -7.0,
        "CDELT1" => -0.02, "CDELT2" => 0.035,
        "PC1_1" => 0.8, "PC1_2" => 0.3, "PC2_1" => -0.2, "PC2_2" => 1.1,
        "CPDIS1" => "TPD", "CPDIS2" => "TPD",
        "DP1.NAXES" => 2, "DP1.AXIS.1" => 1, "DP1.AXIS.2" => 2,
        "DP1.OFFSET.1" => 1.0, "DP1.OFFSET.2" => -2.0,
        "DP1.SCALE.1" => 0.5, "DP1.SCALE.2" => 2.0,
        "DP1.TPD.FWD.0" => 0.01, "DP1.TPD.FWD.1" => 1.02, "DP1.TPD.FWD.2" => -0.03,
        "DP1.TPD.FWD.4" => 0.0007, "DP1.TPD.FWD.5" => -0.0002, "DP1.TPD.FWD.6" => 0.0005,
        "DP2.NAXES" => 2, "DP2.AXIS.1" => 2, "DP2.AXIS.2" => 1,
        "DP2.TPD.FWD.4" => -0.002, "DP2.TPD.FWD.5" => 0.0003,
    )
    prior_tpd = parse_fits_wcs(prior_tpd_header)
    prior_world = pixel_to_world(prior_tpd, 12.3, -3.4)
    @test prior_world ≈ SVector(5.00963112, -7.95957814925) atol=2e-14
    @test world_to_pixel(prior_tpd, prior_world) ≈ SVector(12.3, -3.4) atol=2e-12
    pixel_to_world(prior_tpd, SVector(12.3, -3.4))
    @test allocated_fits_forward(prior_tpd, SVector(12.3, -3.4)) == 0

    sequent_tpd_header = Dict(
        "NAXIS" => 2, "CTYPE1" => "LINEAR", "CTYPE2" => "LINEAR",
        "CRPIX1" => 10.0, "CRPIX2" => 20.0, "CRVAL1" => 5.0, "CRVAL2" => -7.0,
        "CDELT1" => -0.02, "CDELT2" => 0.035,
        "PC1_1" => 0.8, "PC1_2" => 0.3, "PC2_1" => -0.2, "PC2_2" => 1.1,
        "CQDIS1" => "TPD", "CQDIS2" => "TPD",
        "DQ1.NAXES" => 2, "DQ1.AXIS.1" => 1, "DQ1.AXIS.2" => 2,
        "DQ1.OFFSET.1" => 1.0, "DQ1.OFFSET.2" => -2.0,
        "DQ1.SCALE.1" => 0.5, "DQ1.SCALE.2" => 2.0, "DQ1.DOCORR" => 0,
        "DQ1.AUX.1.COEFF.0" => 0.1, "DQ1.AUX.1.COEFF.1" => 1.1, "DQ1.AUX.1.COEFF.2" => -0.2,
        "DQ1.AUX.2.COEFF.0" => -0.3, "DQ1.AUX.2.COEFF.1" => 0.4, "DQ1.AUX.2.COEFF.2" => 0.9,
        "DQ1.TPD.FWD.0" => 0.01, "DQ1.TPD.FWD.1" => 1.02, "DQ1.TPD.FWD.2" => -0.03,
        "DQ1.TPD.FWD.4" => 0.0007, "DQ1.TPD.FWD.5" => -0.0002, "DQ1.TPD.FWD.6" => 0.0005,
        "DQ1.TPD.REV.1" => 1.0,
        "DQ2.NAXES" => 2, "DQ2.AXIS.1" => 2, "DQ2.AXIS.2" => 1,
        "DQ2.TPD.FWD.4" => -0.002, "DQ2.TPD.FWD.5" => 0.0003,
    )
    sequent_tpd = parse_fits_wcs(sequent_tpd_header)
    sequent_world = pixel_to_world(sequent_tpd, 12.3, -3.4)
    @test sequent_world ≈ SVector(4.820512437282, -7.963625782) atol=2e-14
    @test world_to_pixel(sequent_tpd, sequent_world; initial_guess=(12.3, -3.4)) ≈ SVector(12.3, -3.4) atol=2e-12
    pixel_to_world(sequent_tpd, SVector(12.3, -3.4))
    @test allocated_fits_forward(sequent_tpd, SVector(12.3, -3.4)) == 0
    # Prior and sequent TPD keep a third, unrelated axis separable.
    spectral_axis = Dict("NAXIS" => 3, "CTYPE3" => "WAVE", "CRPIX3" => 1.0, "CRVAL3" => 5e-7, "CDELT3" => 1e-10)
    for header in (prior_tpd_header, sequent_tpd_header)
        @test axis_correlation_matrix(parse_fits_wcs(merge(header, spectral_axis))) == SMatrix{3,3,Bool}([1 1 0; 1 1 0; 0 0 1])
    end

    # WCSLIB 8.9 and Astropy 8.0.1/WCSLIB 8.6 references for combined AUX and radial TPD.
    auxiliary_radial_header = Dict(
        "NAXIS" => 2, "NAXIS1" => 100, "NAXIS2" => 100, "CTYPE1" => "LINEAR", "CTYPE2" => "LINEAR",
        "CRPIX1" => 0.0, "CRPIX2" => 0.0, "CRVAL1" => 0.0, "CRVAL2" => 0.0, "CDELT1" => 1.0, "CDELT2" => 1.0,
        "CPDIS1" => "TPD", "CPDIS2" => "TPD", "DP1.NAXES" => 2, "DP1.AXIS.1" => 1, "DP1.AXIS.2" => 2,
        "DP1.AUX.1.COEFF.0" => 1.0, "DP1.AUX.1.COEFF.1" => 2.0, "DP1.AUX.1.COEFF.2" => -0.5,
        "DP1.AUX.2.COEFF.0" => -2.0, "DP1.AUX.2.COEFF.1" => 0.25, "DP1.AUX.2.COEFF.2" => 1.5,
        "DP1.TPD.FWD.3" => 0.1, "DP2.DOCORR" => 0, "DP2.NAXES" => 2, "DP2.AXIS.1" => 1, "DP2.AXIS.2" => 2,
        "DP2.TPD.FWD.2" => 1.0,
    )
    direct_sip_header = Dict(
        "NAXIS" => 2, "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CRPIX1" => 0.0, "CRPIX2" => 0.0,
        "CRVAL1" => 0.0, "CRVAL2" => 0.0, "CDELT1" => 1.0, "CDELT2" => 1.0,
        "CPDIS1" => "SIP", "DP1.NAXES" => 2, "DP1.SIP.FWD.2_0" => 0.001, "DP1.SIP.REV.2_0" => -0.001,
        "CPDIS2" => "SIP", "DP2.NAXES" => 2, "DP2.SIP.FWD.0_2" => 0.002, "DP2.SIP.REV.0_2" => -0.002,
    )
    direct_sip = parse_fits_wcs(direct_sip_header)
    # WCSLIB 8.9 and Astropy 6.1.7/WCSLIB 8.3 references.
    direct_sip_references = (
        (SVector(0.0, 0.0), SVector(0.0, 0.0)),
        (SVector(2.0, 3.0), SVector(2.0031834025537707, 3.0133741951373705)),
        (SVector(10.0, -5.0), SVector(9.9972923518987233, -4.8631269229232190)),
    )
    for (pixel, world) in direct_sip_references
        @test pixel_to_world(direct_sip, pixel) ≈ world atol=5e-14
        iszero(pixel) || @test world_to_pixel(direct_sip, world) ≈ pixel atol=5e-12
    end
    @test_throws "requires DP1.NAXES" parse_fits_wcs(delete!(copy(direct_sip_header), "DP1.NAXES"))
    @test_throws "must be 2 for SIP" parse_fits_wcs(merge(direct_sip_header, Dict("DP1.NAXES" => 1)))
    @test_throws "invalid SIP powers" parse_fits_wcs(merge(direct_sip_header, Dict("DP1.SIP.FWD.9_1" => 1.0)))
    @test_throws "unrecognized SIP parameter" parse_fits_wcs(merge(direct_sip_header, Dict("DP1.DOCORR" => 0)))

    auxiliary_radial = parse_fits_wcs(auxiliary_radial_header)
    auxiliary_radial_references = (
        (SVector(4.0, 3.0), SVector(4.8276472678623428, 3.0)),
        (SVector(-2.0, 5.0), SVector(-1.2566965626340747, 5.0)),
        (SVector(10.0, -4.0), SVector(12.364846718077093, -4.0)),
    )
    for (pixel, world) in auxiliary_radial_references
        @test pixel_to_world(auxiliary_radial, pixel) ≈ world atol=5e-13 rtol=5e-14
        @test world_to_pixel(auxiliary_radial, world) ≈ pixel atol=5e-12 rtol=5e-13
    end

    # WCSLIB 8.9 and Astropy 8.0.1/WCSLIB 8.6 references for a quadratic ZPN.
    quadratic_zpn_header = Dict(
        "NAXIS" => 2, "NAXIS1" => 100, "NAXIS2" => 100, "CTYPE1" => "RA---ZPN", "CTYPE2" => "DEC--ZPN",
        "CRPIX1" => 0.0, "CRPIX2" => 0.0, "CRVAL1" => 10.0, "CRVAL2" => 20.0, "CDELT1" => 1.0, "CDELT2" => 1.0,
        "LONPOLE" => 180.0, "PV2_0" => 0.0, "PV2_1" => 1.0, "PV2_2" => -0.2,
    )
    quadratic_zpn = parse_fits_wcs(quadratic_zpn_header)
    quadratic_zpn_references = (
        (SVector(10.0, 0.0), SVector(21.025845188580881, 19.659355275051063)),
        (SVector(50.0, 0.0), SVector(75.896052659322834, 8.4547122997169151)),
        (SVector(70.0, 0.0), SVector(130.12786649903811, -10.353087255643118)),
    )
    for (pixel, world) in quadratic_zpn_references
        @test pixel_to_world(quadratic_zpn, pixel) ≈ world atol=1e-11 rtol=5e-14
        @test world_to_pixel(quadratic_zpn, world) ≈ pixel atol=5e-11 rtol=5e-13
    end
    @test all(isnan, pixel_to_world(quadratic_zpn, 72.0, 0.0))
    @test_throws DomainError pixel_to_world(quadratic_zpn, 72.0, 0.0; strict=true)
    # WCSLIB defaults ZPN coefficients to zero and rejects these degenerate polynomials.
    bare_zpn_header = filter(pair -> !startswith(pair.first, "PV2_"), quadratic_zpn_header)
    @test_throws "must have a nonzero coefficient" parse_fits_wcs(bare_zpn_header)
    @test_throws "must have a positive linear coefficient" parse_fits_wcs(merge(bare_zpn_header, Dict("PV2_3" => 1e-3)))
    skipped_tpd = parse_fits_wcs(merge(sequent_tpd_header, Dict("CQERR1" => 0.01, "CQERR2" => 0.01)); minerr=0.02)
    @test pixel_to_world(skipped_tpd, 12.3, -3.4) ≈ SVector(5.1036, -7.917) atol=2e-15

    reverse_tpd_header = Dict(
        "NAXIS" => 2, "CTYPE1" => "LINEAR", "CTYPE2" => "LINEAR", "CRPIX1" => 0.0, "CRPIX2" => 0.0,
        "CRVAL1" => 0.0, "CRVAL2" => 0.0, "CDELT1" => 1.0, "CDELT2" => 1.0,
        "CQDIS1" => "TPD", "CQDIS2" => "TPD", "DQ1.NAXES" => 2, "DQ2.NAXES" => 2,
        "DQ1.TPD.FWD.1" => 1.0, "DQ1.TPD.FWD.4" => 0.1, "DQ1.TPD.REV.4" => -0.08,
        "DQ2.TPD.FWD.2" => 1.0,
    )
    reverse_tpd = parse_fits_wcs(reverse_tpd_header)
    reverse_world = pixel_to_world(reverse_tpd, 2.0, 3.0)
    reverse_estimate, used_reverse = GWCS._inverse_initialize(reverse_tpd.forward, reverse_world)
    @test used_reverse
    @test reverse_estimate ≈ SVector(2.8512, 6.0)
    @test maximum(abs, pixel_to_world(reverse_tpd, reverse_estimate) - reverse_world) > 1
    @test world_to_pixel(reverse_tpd, reverse_world) ≈ SVector(2.0, 3.0) atol=1e-12
    @test world_to_pixel(reverse_tpd, reverse_world; initial_guess=(2.0, 3.0), maxiter=1) ≈ SVector(2.0, 3.0) atol=1e-12

    high_order_tpd_header = Dict(
        "NAXIS" => 2, "CTYPE1" => "LINEAR", "CTYPE2" => "LINEAR",
        "CRPIX1" => 0.0, "CRPIX2" => 0.0, "CRVAL1" => 0.0, "CRVAL2" => 0.0,
        "CDELT1" => 1.0, "CDELT2" => 1.0,
        "CPDIS1" => "TPD", "DP1.NAXES" => 1, "DP1.AXIS.1" => 1,
        "DP1.TPD.FWD.49" => 1e-9, "DP1.TPD.FWD.59" => 10.0,
    )
    @test pixel_to_world(parse_fits_wcs(high_order_tpd_header), 2.0, 3.0) ≈ SVector(2.000000512, 3.0) atol=2e-15
    radial_tpd_header = merge(high_order_tpd_header, Dict("DP1.NAXES" => 2, "DP1.AXIS.2" => 2, "DP1.TPD.FWD.49" => 0.0, "DP1.TPD.FWD.59" => 1e-9))
    @test pixel_to_world(parse_fits_wcs(radial_tpd_header), 3.0, 4.0) ≈ SVector(3.001953125, 4.0) atol=2e-15

    spectral_tpd_header = Dict(
        "NAXIS" => 1, "CTYPE1" => "VOPT-F2W", "CUNIT1" => "m/s",
        "CRPIX1" => 1.0, "CRVAL1" => 1000.0, "CDELT1" => 100.0, "RESTWAV" => 0.21,
        "CQDIS1" => "TPD", "DQ1.NAXES" => 1, "DQ1.AXIS.1" => 1, "DQ1.TPD.FWD.1" => 0.1,
    )
    @test pixel_to_world(parse_fits_wcs(spectral_tpd_header), 2.5) ≈ SVector(1165.000090777874) atol=5e-8

    legacy_header(path) = Dict{String,Any}(begin
        result = Pair{String,Any}[]
        for line in eachline(path)
            startswith(line, '#') && continue
            length(line) < 10 && continue
            equal = findfirst(==('='), line)
            equal === nothing && continue
            keyword = strip(line[1:8])
            raw = strip(first(split(line[(equal + 1):end], '/'; limit=2)))
            value = startswith(raw, "'") ? strip(raw, ['\'', ' ']) : something(tryparse(Int, raw), something(tryparse(Float64, raw), raw))
            push!(result, keyword => value)
        end
        result
    end)
    legacy_references = (
        "DSS" => SVector(81.938890331523893, -12.765680207379495),
        "TNX" => SVector(266.73116925327901, -30.13403790309048),
        "ZPX" => SVector(321.03938306885675, 37.191778590405455),
    )
    for (kind, expected) in legacy_references
        legacy = parse_fits_wcs(legacy_header(joinpath(@__DIR__, "fixtures", "wcslib", "$kind.keyrec")))
        world = pixel_to_world(legacy, 200.0, 200.0)
        @test world ≈ expected atol=2e-12
        @test world_to_pixel(legacy, world) ≈ SVector(200.0, 200.0) atol=2e-7
        @test @inferred(pixel_to_world(legacy, SVector(200.0, 200.0))) ≈ expected atol=2e-12
        pixel_to_world(legacy, SVector(200.0, 200.0))
        @test allocated_fits_forward(legacy, SVector(200.0, 200.0)) == 0
    end
    dss_with_cd = legacy_header(joinpath(@__DIR__, "fixtures", "wcslib", "DSS.keyrec"))
    dss_with_cd["CD1_1"] = 1.0
    dss_with_cd["CD1_2"] = dss_with_cd["CD2_1"] = 0.0
    dss_with_cd["CD2_2"] = 1.0
    @test pixel_to_world(parse_fits_wcs(dss_with_cd), 200.0, 200.0) ≈ last(first(legacy_references)) atol=2e-12
    @test_throws "requires both WAT1_n and WAT2_n" parse_fits_wcs(merge(tan_header, Dict("CTYPE1" => "RA---TNX", "CTYPE2" => "DEC--TNX", "WAT1_001" => "wtype=tnx")))
    # Like WCSLIB, WAT strings that do not encode TNX/ZPX (common in IRAF spectra) are ignored.
    linear_wat_header = Dict("NAXIS" => 1, "CTYPE1" => "LINEAR", "CRPIX1" => 1.0, "CRVAL1" => 4000.0, "CDELT1" => 1.5,
        "WAT0_001" => "system=equispec", "WAT1_001" => "wtype=linear label=Wavelength units=angstroms")
    @test pixel_to_world(parse_fits_wcs(linear_wat_header), 3.0) == SVector(4003.0)
    tan_wat_header = merge(tan_header, Dict("WAT0_001" => "system=image", "WAT1_001" => "wtype=tan axtype=ra", "WAT2_001" => "wtype=tan axtype=dec"))
    @test pixel_to_world(parse_fits_wcs(tan_wat_header), 400.0, 300.0) == pixel_to_world(tan, 400.0, 300.0)
    @test_throws "unsupported WAT1 wtype" parse_fits_wcs(merge(tan_wat_header, Dict("CTYPE1" => "RA---TNX", "CTYPE2" => "DEC--TNX")))

    # Like astropy, -TAN-SIP with SCAMP PV terms uses SIP and ignores the celestial PV keywords,
    # while plain -TAN with SCAMP PV terms is TPV and ignores SIP.
    scamp_header = merge(sip_header, Dict("PV1_1" => 1.0, "PV1_5" => 1e-6, "PV2_1" => 1.0))
    original_scamp_header = copy(scamp_header)
    scamp = @test_logs (:warn, r"Ignoring SCAMP PVi_m") parse_fits_wcs(scamp_header)
    @test pixel_to_world(scamp, 500.0, 500.0) ≈ pixel_to_world(sip, 500.0, 500.0) atol=1e-13
    @test scamp_header == original_scamp_header
    pre2012_scamp = merge(scamp_header, Dict("CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN"))
    pre2012_without_sip = filter(pair -> !startswith(pair.first, "A_") && !startswith(pair.first, "B_"), pre2012_scamp)
    @test pixel_to_world(parse_fits_wcs(pre2012_scamp), 500.0, 500.0) ≈ pixel_to_world(parse_fits_wcs(pre2012_without_sip), 500.0, 500.0) atol=1e-13
    @test pixel_to_world(parse_fits_wcs(pre2012_scamp), 500.0, 500.0) != pixel_to_world(sip, 500.0, 500.0)

    projection_references = (
        "TAN" => (29.89839490759596, 10.199983486206756), "SIN" => (29.89839413334517, 10.19998500905114),
        "STG" => (29.898394520473516, 10.199984247623151), "ARC" => (29.89839439143152, 10.199984501430938),
        "ZEA" => (29.898394326910267, 10.199984628335333), "AZP" => (29.89839490759596, 10.199983486206756),
        "SZP" => (29.898394907596106, 10.199983486206472), "CAR" => (29.89839480407761, 10.19998460309562),
        "CEA" => (29.898410445812885, 10.14998477681374), "CYP" => (29.898394804141162, 10.199984400017723),
        "MER" => (29.898394804204713, 10.19998419694032), "SFL" => (29.898394185060283, 10.199984602908012),
        "PAR" => (29.898397122608117, 10.19097057437874), "MOL" => (29.88715158733106, 10.180044382679622),
        "PCO" => (29.89839418515686, 10.199984298289833), "AIT" => (29.89839433656448, 10.199984685454476),
        "ZPN" => (29.89839439143152, 10.199984501430938), "COP" => (29.898394804151437, 10.19998379183283),
        "COD" => (29.898394805341276, 10.199984603072604), "COE" => (29.898394805937766, 10.199985009762853),
        "COO" => (29.898394804745834, 10.199984197096427), "BON" => (29.898394184440406, 10.199935629296364),
        "TSC" => (29.87061068884833, 10.254620638891618), "CSC" => (29.90594452860057, 10.185134498675133),
        "QSC" => (29.912650673686358, 10.187270730775742), "HPX" => (29.8984042630899, 10.169750125560459),
        "XPH" => (29.887966865612718, 10.190967211678425),
    )
    for (code, expected) in projection_references
        header = Dict(
            "NAXIS" => 2, "CTYPE1" => "RA---$code", "CTYPE2" => "DEC--$code",
            "CRPIX1" => 10.0, "CRPIX2" => 20.0, "CRVAL1" => 30.0, "CRVAL2" => 10.0,
            "CDELT1" => -0.1, "CDELT2" => 0.1,
        )
        code == "CEA" && (header["PV2_1"] = 0.75)
        code == "ZPN" && (header["PV2_1"] = 1.0)
        code in ("COP", "COD", "COE", "COO", "BON") && (header["PV2_1"] = 45.0)
        world = pixel_to_world(parse_fits_wcs(header), 11.0, 22.0)
        tolerance = code == "CSC" ? 1e-8 : 5e-12
        @test world ≈ SVector(expected) atol=tolerance
    end

    # WCSLIB 8.10 references; with CRPIX = 0 and CDELT = 1 the intermediate coordinates equal the pixel coordinates.
    unit_projection_header(code, crval, extra=Dict()) = merge(Dict("NAXIS" => 2, "CTYPE1" => "RA---$code", "CTYPE2" => "DEC--$code",
        "CRPIX1" => 0.0, "CRPIX2" => 0.0, "CRVAL1" => crval[1], "CRVAL2" => crval[2], "CDELT1" => 1.0, "CDELT2" => 1.0), extra)
    southern_references = (
        "COP" => (40.17098347253864, -19.978900355365333), "COD" => (40.19477387030446, -19.43260045551106),
        "COE" => (40.22339435934347, -19.347847108380996), "COO" => (40.16487877374618, -19.544565363066575),
        "BON" => (40.5841787531637, -18.958216377056456),
    )
    for (code, expected) in southern_references
        parameters = code == "BON" ? Dict("PV2_1" => -45.0) : Dict("PV2_1" => -45.0, "PV2_2" => 10.0)
        southern = parse_fits_wcs(unit_projection_header(code, (30.0, -40.0), parameters))
        @test pixel_to_world(southern, 10.0, 20.0) ≈ SVector(expected) atol=1e-12
        @test world_to_pixel(southern, expected...) ≈ SVector(10.0, 20.0) atol=1e-10
    end
    # WCSLIB requires PV2_1 for the conic and Bonne projections rather than defaulting it.
    for code in ("COP", "COD", "COE", "COO", "BON")
        @test_throws "FITS $code projection requires PV2_1, as in WCSLIB" parse_fits_wcs(unit_projection_header(code, (30.0, 40.0)))
    end
    # Points beyond the projection boundary have no sky coordinates (WCSLIB status 1), not wrapped duplicates.
    for (code, parameters, pixel) in (("CAR", Dict(), (-189.0, 0.0)), ("MOL", Dict(), (170.0, 60.0)), ("SFL", Dict(), (170.0, 60.0)),
            ("HPX", Dict(), (-170.0, 80.0)), ("XPH", Dict(), (99.0, -138.0)), ("COD", Dict("PV2_1" => 45.0, "PV2_2" => 25.0), (195.0, -156.0)))
        @test all(isnan, pixel_to_world(parse_fits_wcs(unit_projection_header(code, (0.0, 0.0), parameters)), pixel...))
        @test_throws DomainError pixel_to_world(parse_fits_wcs(unit_projection_header(code, (0.0, 0.0), parameters)), pixel...; strict=true)
    end
    # Polyconic intermediate y exceeds 90 degrees away from the central meridian.
    @test pixel_to_world(parse_fits_wcs(unit_projection_header("PCO", (0.0, 0.0))), 67.0, -132.0) ≈ SVector(172.10045640795752, -39.05258488360962) atol=1e-11
end

@testset "Additional shared spherical projections" begin
    projections = (
        PlateCarreeProjection(), CylindricalEqualArea(), CylindricalPerspective(), MercatorProjection(),
        SansonFlamsteedProjection(), ParabolicProjection(), MollweideProjection(), PolyconicProjection(), HammerAitoffProjection(),
        ZenithalPolynomial(), ConicPerspective(), ConicEquidistant(), ConicEqualArea(), ConicOrthomorphic(), BonneProjection(),
        TangentialSphericalCube(), QuadrilateralizedSphericalCube(), HEALPixProjection(), HEALPixPolarProjection(),
    )
    native = SVector(20.0, 30.0)
    for projection in projections
        forward = NativeToIntermediate(projection)
        @test inverse_transform(forward)(forward(native)) ≈ native atol=5e-11
        @test eltype(forward(SVector(20.0f0, 30.0f0))) === Float32
        @test transform_dependencies(forward) == @SMatrix [true true; true true]
        @test size(transform_jacobian(forward, native)) == (2, 2)
    end
    cobe = NativeToIntermediate(COBESphericalCube())
    @test inverse_transform(cobe)(cobe(native)) ≈ native atol=3e-3
end

# References are WCSLIB 8.10 (where it accepts the header) and Astropy 8.0.1, which applies WCSLIB's wcsfix repairs.
@testset "Legacy FITS conventions" begin
    celestial(code; crval=(30.0, 40.0), extra...) = merge(Dict{String,Any}("NAXIS" => 2, "CTYPE1" => "RA---$code", "CTYPE2" => "DEC--$code",
        "CRPIX1" => 50.0, "CRPIX2" => 50.0, "CRVAL1" => crval[1], "CRVAL2" => crval[2], "CDELT1" => -0.1, "CDELT2" => 0.1), Dict{String,Any}(String(k) => v for (k, v) in extra))
    spectral(ctype, crval, cdelt; extra...) = merge(Dict{String,Any}("NAXIS" => 1, "CTYPE1" => ctype, "CRPIX1" => 32.0, "CRVAL1" => crval, "CDELT1" => cdelt,
        "RESTFRQ" => 1.420405751e9), Dict{String,Any}(String(k) => v for (k, v) in extra))
    repaired(header, pattern; kwargs...) = @test_logs (:warn, pattern) parse_fits_wcs(header; kwargs...)

    ncp = celestial("NCP")
    @test pixel_to_world(repaired(ncp, r"NCP projection read as SIN"), 60.0, 43.0) ≈ SVector(28.708058345297673, 39.27700409069822) atol=1e-12
    @test_throws "NCP projection is undefined for a zero reference latitude" parse_fits_wcs(celestial("NCP"; crval=(30.0, 0.0)))
    @test pixel_to_world(repaired(celestial("GLS"), r"GLS projection read as SFL"), 60.0, 43.0) ≈ SVector(28.70774355295724, 39.3) atol=1e-12
    @test pixel_to_world(repaired(merge(ncp, Dict("CTYPE1" => "RA", "CTYPE2" => "DEC")), r"no projection code"), 60.0, 43.0) ≈ SVector(29.0, 39.3) atol=1e-12

    # AIPS spectral types: the suffix (or VELREF) gives SPECSYS; VELO is optical unless VELREF/256 == 1.
    frame(wcs) = output_frame(wcs).reference_frame.spectral.specsys
    @test frame(repaired(spectral("FREQ-LSR", 1.42e9, -2e4), r"\"FREQ-LSR\" read as \"FREQ\" with SPECSYS=\"LSRK\"")) == "LSRK"
    velref_freq = repaired(spectral("FREQ-OBS", 1.42e9, -2e4; VELREF=2), r"VELREF=2")
    @test frame(velref_freq) == "BARYCENT"
    @test pixel_to_world(velref_freq, 40.0) ≈ SVector(1.41984e9)
    @test pixel_to_world(repaired(spectral("FELO-HEL", 1.5e4, 4e3), r"read as \"VOPT-F2W\""), 40.0) ≈ SVector(47003.415890040655) rtol=1e-12
    optical = repaired(spectral("VELO-LSR", 1.5e4, 4e3), r"read as \"VOPT\"")
    radio = repaired(spectral("VELO-LSR", 1.5e4, 4e3; VELREF=257), r"read as \"VRAD\"")
    @test pixel_to_world(optical, 40.0) == pixel_to_world(radio, 40.0) == SVector(47000.0)
    @test world_axis_physical_types(optical) == ("spect.dopplerVeloc.opt",)
    @test world_axis_physical_types(radio) == ("spect.dopplerVeloc.radio",)
    @test frame(parse_fits_wcs(spectral("VELO-LSR", 1.5e4, 4e3; SPECSYS="SOURCE"); legacy=:quiet)) == "SOURCE"
    @test_throws "invalid AIPS reference frame code 9" parse_fits_wcs(spectral("FREQ-LSR", 1.42e9, -2e4; VELREF=9))
    @test_throws "unsupported nonlinear algorithm" parse_fits_wcs(spectral("FREQ-XYZ", 1.42e9, -2e4; VELREF=1))

    # Keyword spellings from early drafts of Papers I-III.
    reference = parse_fits_wcs(merge(ncp, Dict("CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "PC1_1" => 0.9, "PC1_2" => 0.2, "PC2_1" => -0.1, "PC2_2" => 1.1)))
    @test pixel_to_world(reference, 60.0, 43.0) ≈ SVector(29.0204675322597, 39.125967753957475) atol=1e-12
    for padded in (Dict("PC001001" => 0.9, "PC001002" => 0.2, "PC002001" => -0.1, "PC002002" => 1.1), Dict("PC01_01" => 0.9, "PC01_02" => 0.2, "PC02_01" => -0.1, "PC02_02" => 1.1))
        @test pixel_to_world(repaired(merge(celestial("TAN"), padded), r"read as PC1_1"), 60.0, 43.0) == pixel_to_world(reference, 60.0, 43.0)
    end
    padded_cd = filter(pair -> !startswith(pair.first, "CDELT"), merge(celestial("TAN"), Dict("CD001001" => -0.09, "CD001002" => 0.02, "CD002001" => 0.01, "CD002002" => 0.11)))
    @test pixel_to_world(repaired(padded_cd, r"CD001001 read as CD1_1"), 60.0, 43.0) ≈ SVector(28.65581489581939, 39.322300146818115) atol=1e-12
    @test_throws "PC001001 and PC1_1 give conflicting values" parse_fits_wcs(merge(celestial("TAN"), Dict("PC001001" => 0.9, "PC1_1" => 1.0)))
    @test pixel_to_world(repaired(celestial("SIN"; PROJP1=0.1, PROJP2=-0.2), r"PROJP1 read as PV2_1"), 60.0, 43.0) ≈ SVector(28.70604087434814, 39.29531321986866) atol=1e-12
    epoch = output_frame(repaired(celestial("TAN"; EPOCH=1950.0), r"EPOCH read as EQUINOX")).reference_frame
    @test (epoch.radesys, epoch.equinox) == ("FK4", 1950.0)
    @test output_frame(parse_fits_wcs(celestial("TAN"; EPOCH=1950.0, EQUINOX=2000.0))).reference_frame.equinox == 2000.0
    @test output_frame(repaired(celestial("TAN"; RADECSYS="FK5"), r"RADECSYS read as RADESYS")).reference_frame.radesys == "FK5"
    @test output_frame(repaired(spectral("VRAD", 1.5e4, 4e3; VSOURCE=1.5e4), r"VSOURCE read as ZSOURCE")).reference_frame.spectral.zsource ≈ 5.0035866073772084e-5 rtol=1e-12
    time_header = Dict{String,Any}("NAXIS" => 1, "CTYPE1" => "TIME", "CRPIX1" => 1.0, "CRVAL1" => 0.0, "CDELT1" => 1.0)
    epoch_of(header) = output_frame(parse_fits_wcs(header; legacy=:quiet)).reference_frame.time.reference_epoch
    @test output_frame(repaired(merge(time_header, Dict("MJD-REF" => 59000.5)), r"MJD-REF read as MJDREF")).reference_frame.time.reference_epoch == 59000.5
    @test epoch_of(merge(time_header, Dict("JDREF" => 2459000.5))) == 59000.0
    @test epoch_of(merge(time_header, Dict("JDREFI" => 2459000, "JDREFF" => 0.75))) == 59000.25

    # wcsfix repairs that astropy applies by default.
    cube = Dict{String,Any}("NAXIS" => 3, "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CTYPE3" => "FREQ", "CRPIX1" => 25.0, "CRPIX2" => 25.0, "CRPIX3" => 1.0,
        "CRVAL1" => 30.0, "CRVAL2" => 40.0, "CRVAL3" => 1.4e9, "CD1_1" => -0.01, "CD2_2" => 0.01, "CDELT3" => 1e6)
    @test pixel_to_world(repaired(cube, r"CD3_3 set to 1"), 30.0, 20.0, 5.0) ≈ SVector(29.934777423125915, 39.94998173882473, 1.400000004e9) atol=1e-12
    @test_throws "matrix from CDi_j is singular" parse_fits_wcs(merge(celestial("TAN"), Dict("CD1_1" => 1.0, "CD1_2" => 1.0, "CD2_1" => 1.0, "CD2_2" => 1.0)))
    unsuffixed_sip = celestial("TAN"; A_ORDER=2, B_ORDER=2, A_2_0=1e-5, B_0_2=1e-5)
    @test pixel_to_world(repaired(unsuffixed_sip, r"SIP coefficients applied although CTYPE lacks the -SIP suffix"), 60.0, 43.0) ≈ SVector(28.707928277445472, 39.2929427522218) atol=1e-12
    @test pixel_to_world(repaired(celestial("TAN"; CROTA1=20.0), r"CROTA1 on the longitude axis is ignored"), 60.0, 43.0) == pixel_to_world(parse_fits_wcs(celestial("TAN")), 60.0, 43.0)

    # Units are case-sensitive; only WCSLIB's safe translations are applied, as in astropy.
    @test pixel_to_world(repaired(spectral("WAVE", 6563.0, 0.5; CUNIT1="Angstroms"), r"\"Angstroms\" read as \"Angstrom\""), 40.0) ≈ SVector(6.567e-7) rtol=1e-14
    @test pixel_to_world(repaired(spectral("FREQ", 1420.0, -0.02; CUNIT1="MHZ"), r"\"MHZ\" read as \"MHz\""), 40.0) ≈ SVector(1.41984e9)
    @test world_axis_units(parse_fits_wcs(spectral("FREQ", 1420.0, -0.02; CUNIT1="MHZ"); preserve_units=true, legacy=:quiet)) == ("MHz",)
    @test_throws "unsupported spectral CUNIT1 (written \"KM/S\")=\"km/S\"" parse_fits_wcs(spectral("VRAD", 15.0, 4.0; CUNIT1="KM/S"); legacy=:quiet)
    @test_throws "unsupported spectral CUNIT1=\"mhz\"" parse_fits_wcs(spectral("FREQ", 1420.0, -0.02; CUNIT1="mhz"))
    @test_throws "unsupported spectral CUNIT1=\"micron\"" parse_fits_wcs(spectral("WAVE", 0.6563, 1e-5; CUNIT1="micron"))
    @test_throws "ambiguous time CUNIT1=\"S\"" parse_fits_wcs(merge(time_header, Dict("CUNIT1" => "S")))
    @test pixel_to_world(repaired(merge(time_header, Dict("CUNIT1" => "SEC")), r"\"SEC\" read as \"s\""), 3.0) == SVector(2.0)

    # Repair modes and alternate keys.
    @test_logs parse_fits_wcs(ncp; legacy=:quiet)
    @test_throws "uses legacy conventions" parse_fits_wcs(ncp; legacy=:error)
    @test_throws "legacy must be :warn, :quiet, or :error" parse_fits_wcs(ncp; legacy=:strict)
    @test parse_fits_wcs(celestial("TAN"); legacy=:error) isa WCS
    alternate = Dict{String,Any}("NAXIS" => 1, "WCSAXESA" => 1, "CTYPE1A" => "VELO-HEL", "CRPIX1A" => 32.0, "CRVAL1A" => 1.5e4, "CDELT1A" => 4e3, "VELREFA" => 257, "EPOCHA" => 1950.0)
    alternate_wcs = repaired(alternate, r"key 'A'"; key='A')
    @test world_axis_physical_types(alternate_wcs) == ("spect.dopplerVeloc.radio",)
    # VELREF's frame code (257 mod 256 = 1, LSRK) takes precedence over the -HEL suffix, as in WCSLIB.
    @test frame(alternate_wcs) == "LSRK"
    @test_throws "FITS WCS key 'A' is malformed" find_all_fits_wcs(alternate; legacy=:error)
    legacy_cards = reduce(vcat, FITSFiles.Card(key, value) for (key, value) in ("CTYPE1" => "RA---NCP", "CTYPE2" => "DEC--NCP", "CRPIX1" => 50.0, "CRPIX2" => 50.0,
        "CRVAL1" => 30.0, "CRVAL2" => 40.0, "CDELT1" => -0.1, "CDELT2" => 0.1))
    @test_throws "uses legacy conventions" parse_fits_wcs(legacy_cards; legacy=:error)
    @test pixel_to_world(parse_fits_wcs(legacy_cards; legacy=:quiet), 60.0, 43.0) ≈ SVector(28.708058345297673, 39.27700409069822) atol=1e-12
end
