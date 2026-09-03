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
    @test_throws "must not be mixed" parse_fits_wcs(merge(pc_header, Dict("CD1_1" => 1.0)))

    crota_header = merge(linear_header, Dict("CDELT1" => 1.0, "CDELT2" => 1.0, "CROTA2" => 90.0))
    @test pixel_to_world(parse_fits_wcs(crota_header), 11.0, 20.0) ≈ SVector(100.0, 201.0) atol=1e-14

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

    time_header = Dict(
        "NAXIS" => 1, "CTYPE1" => "TIME", "TIMEUNIT" => "d", "TIMESYS" => "UTC", "TREFPOS" => "TOPOCENTER",
        "MJDREFI" => 59000, "MJDREFF" => 0.25, "OBSGEO-X" => 1.0, "OBSGEO-Y" => 2.0, "OBSGEO-Z" => 3.0,
        "CRPIX1" => 1.0, "CRVAL1" => 1.0, "CDELT1" => 0.5,
    )
    time = parse_fits_wcs(time_header)
    time_metadata = output_frame(time).reference_frame.time
    @test pixel_to_world(time, 1.0) == SVector(86400.0)
    @test pixel_to_world(parse_fits_wcs(time_header; preserve_units=true), 2.0) == SVector(1.5)
    @test world_axis_units(time) == ("s",)
    @test world_axis_units(parse_fits_wcs(time_header; preserve_units=true)) == ("d",)
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

    fitsio_header = FITSIO.FITSHeader(["NAXIS", "CTYPE1", "CRVAL1"], Any[1, "X", 2.0], ["", "", ""])
    @test pixel_to_world(parse_fits_wcs(fitsio_header), 0.0) == SVector(2.0)
    fitsfiles_cards = FITSFiles.Card[FITSFiles.Card("NAXIS", 1), FITSFiles.Card("CTYPE1", "X"), FITSFiles.Card("CRVAL1", 3.0)]
    @test pixel_to_world(parse_fits_wcs(fitsfiles_cards), 0.0) == SVector(3.0)

    lookup_cards = FITSFiles.Card[
        FITSFiles.Card("CTYPE1", "RA---TAN"), FITSFiles.Card("CTYPE2", "DEC--TAN"), FITSFiles.Card("CRPIX1", 2.0), FITSFiles.Card("CRPIX2", 2.0),
        FITSFiles.Card("CRVAL1", 30.0), FITSFiles.Card("CRVAL2", 10.0), FITSFiles.Card("CDELT1", -0.1), FITSFiles.Card("CDELT2", 0.1),
        FITSFiles.Card("D2IMDIS1", "LOOKUP"), FITSFiles.Card("D2IMERR1", 0.1), FITSFiles.Card("D2IM1", "AXIS.1: 1"),
    ]
    lookup_primary = FITSFiles.HDU(zeros(4, 4), lookup_cards)
    lookup_extension = FITSFiles.HDU(ones(1, 1), FITSFiles.Card[FITSFiles.Card("EXTNAME", "D2IMARR"), FITSFiles.Card("EXTVER", 1), FITSFiles.Card("CRPIX1", 0.0), FITSFiles.Card("CRPIX2", 0.0)])
    fitsfiles_lookup = parse_fits_wcs([lookup_primary, lookup_extension], 1)
    fitsfiles_plain = parse_fits_wcs(merge(tan_header, Dict("CRPIX1" => 2.0, "CRPIX2" => 2.0, "CRVAL1" => 30.0, "CRVAL2" => 10.0, "CDELT1" => -0.1, "CDELT2" => 0.1)))
    @test pixel_to_world(fitsfiles_lookup, 2.0, 2.0) ≈ pixel_to_world(fitsfiles_plain, 3.0, 2.0) atol=1e-12
    @test pixel_to_world(parse_fits_wcs([lookup_primary, lookup_extension], 1; minerr=0.2), 2.0, 2.0) ≈ pixel_to_world(fitsfiles_plain, 2.0, 2.0) atol=1e-12
    malformed_cards = FITSFiles.Card[FITSFiles.Card("CTYPE1", "X"), FITSFiles.Card("WCSAXESA", 1), FITSFiles.Card("CTYPE1A", "DEC--TAN")]
    @test_throws "FITS WCS key 'A' is malformed" find_all_fits_wcs([FITSFiles.HDU(zeros(2), malformed_cards)], 1)

    tab_primary = FITSFiles.HDU(zeros(3), FITSFiles.Card[
        FITSFiles.Card("CTYPE1", "FREQ-TAB"), FITSFiles.Card("CUNIT1", "GHz"), FITSFiles.Card("CRPIX1", 1.0),
        FITSFiles.Card("CRVAL1", 1.0), FITSFiles.Card("CDELT1", 1.0), FITSFiles.Card("PS1_0", "WCS-TABLE"), FITSFiles.Card("PS1_1", "coord"),
    ])
    tab_extension = FITSFiles.HDU((COORD=[10.0, 20.0, 40.0],), FITSFiles.Card[
        FITSFiles.Card("EXTNAME", "WCS-TABLE"), FITSFiles.Card("EXTVER", 1), FITSFiles.Card("EXTLEVEL", 1),
    ])
    fitsfiles_tab = parse_fits_wcs([tab_primary, tab_extension], 1)
    @test pixel_to_world(fitsfiles_tab, 2.5) == SVector(3.0e10)
    @test world_to_pixel(fitsfiles_tab, 3.0e10) == SVector(2.5)

    mktempdir() do directory
        path = joinpath(directory, "lookup.fits")
        primary_header = FITSIO.FITSHeader(["CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2", "D2IMDIS1", "D2IM1"],
            Any["RA---TAN", "DEC--TAN", 2.0, 2.0, 30.0, 10.0, -0.1, 0.1, "LOOKUP", "AXIS.1: 1"], fill("", 10))
        table_header = FITSIO.FITSHeader(["CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2"], Any[0.0, 0.0, 0.0, 0.0, 1.0, 1.0], fill("", 6))
        FITSIO.FITS(path, "w") do file
            FITSIO.write(file, zeros(4, 4); header=primary_header)
            FITSIO.write(file, ones(1, 1); header=table_header, name="D2IMARR", ver=1)
        end
        FITSIO.FITS(path) do file
            @test pixel_to_world(parse_fits_wcs(file[1]), 2.0, 2.0) ≈ pixel_to_world(fitsfiles_plain, 3.0, 2.0) atol=1e-12
        end

        tab_path = joinpath(directory, "tab.fits")
        tab_image_header = FITSIO.FITSHeader(
            ["CTYPE1", "CUNIT1", "CRPIX1", "CRVAL1", "CDELT1", "PS1_0", "PS1_1", "WCSAXESA", "CTYPE1A"],
            Any["FREQ-TAB", "GHz", 1.0, 1.0, 1.0, "WCS-TABLE", "COORD", 1, "DEC--TAN"], fill("", 9),
        )
        tab_table_header = FITSIO.FITSHeader(["EXTLEVEL"], Any[1], [""])
        FITSIO.FITS(tab_path, "w") do file
            FITSIO.write(file, zeros(3); header=tab_image_header)
            FITSIO.write(file, Dict("COORD" => reshape([10.0, 20.0, 40.0], 3, 1)); header=tab_table_header, name="WCS-TABLE", ver=1)
        end
        FITSIO.FITS(tab_path) do file
            fitsio_tab = parse_fits_wcs(file[1])
            @test pixel_to_world(fitsio_tab, 2.5) == SVector(3.0e10)
            @test world_to_pixel(fitsio_tab, 3.0e10) == SVector(2.5)
            @test_throws "FITS WCS key 'A' is malformed" find_all_fits_wcs(file[1])
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
    sip_world = pixel_to_world(sip, 400.0, 300.0)
    @test sip_world ≈ SVector(83.8533100335119, -5.449965264308421) atol=2e-11
    @test world_to_pixel(sip, sip_world; initial_guess=SVector(400.0, 300.0)) ≈ SVector(400.0, 300.0) atol=1e-8

    sip_inverse_header = merge(sip_header, Dict(
        "AP_ORDER" => 2, "AP_2_0" => -5e-6, "AP_0_2" => -2e-6,
        "BP_ORDER" => 2, "BP_2_0" => -1e-6, "BP_1_1" => -3e-6,
    ))
    sip_inverse = parse_fits_wcs(sip_inverse_header)
    @test has_inverse(sip_inverse.forward)
    @test world_to_pixel(sip_inverse, pixel_to_world(sip_inverse, 510.0, 510.0)) ≈ SVector(510.0, 510.0) atol=1e-6

    sip_cube_header = merge(sip_header, Dict("NAXIS" => 3, "CTYPE3" => "FREQ", "CUNIT3" => "Hz", "CRPIX3" => 1.0, "CRVAL3" => 10.0, "CDELT3" => 2.0))
    @test pixel_to_world(parse_fits_wcs(sip_cube_header), 512.0, 512.0, 2.0) ≈ SVector(83.8221, -5.3911, 12.0) atol=1e-11

    lookup_table = BilinearLookupTable(ones(1, 1); crpix=(0, 0))
    lookup_header = merge(tan_header, Dict("D2IMDIS1" => "LOOKUP", "D2IM1.AXIS.1" => 1))
    lookup_data = (detector_to_image=(lookup_table, nothing), cpdis=(nothing, nothing))
    lookup_wcs = GWCS._parse_fits_wcs(GWCS._fits_header(lookup_header); lookup_data)
    @test pixel_to_world(lookup_wcs, 512.0, 512.0) ≈ pixel_to_world(tan, 513.0, 512.0) atol=1e-12
    pixel_to_world(lookup_wcs, SVector(512.0, 512.0))
    @test allocated_fits_forward(lookup_wcs, SVector(512.0, 512.0)) == 0
    bilinear = BilinearLookupTable([0.0 1.0; 2.0 3.0])
    @test bilinear(SVector(0.5, 0.5)) == SVector(1.5)
    @test_throws "CDELT values must be nonzero" BilinearLookupTable(ones(2, 2); cdelt=(1, 0))

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
    @test_throws DomainError pixel_to_world(tab_wcs, 0.49)
    @test pixel_to_world(GWCS._parse_fits_wcs(GWCS._fits_header(tab_header); preserve_units=true, tabular_data=(tabular,)), 2.5) == SVector(30.0)
    pixel_to_world(tab_wcs, SVector(2.5))
    world_to_pixel(tab_wcs, SVector(3.0e10))
    @test allocated_fits_forward(tab_wcs, SVector(2.5)) == 0
    @test allocated_fits_inverse(tab_wcs, SVector(3.0e10)) == 0
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
        "PV1_1" => 1.0, "PV1_4" => 0.01, "PV2_2" => 1.0, "PV2_6" => -0.02,
    ))
    tpv = parse_fits_wcs(tpv_header)
    tpv_world = pixel_to_world(tpv, 513.0, 514.0)
    @test tpv_world ≈ SVector(83.82181819578481, -5.390827775535175) atol=2e-11
    @test world_to_pixel(tpv, tpv_world; initial_guess=(513.0, 514.0)) ≈ SVector(513.0, 514.0) atol=1e-8

    scamp_header = merge(sip_header, Dict("CTYPE1" => "RA---TAN-SIP", "CTYPE2" => "DEC--TAN-SIP", "PV1_1" => 1.0, "PV1_5" => 1e-6, "PV2_2" => 1.0))
    original_scamp_header = copy(scamp_header)
    scamp = parse_fits_wcs(scamp_header)
    scamp_without_sip = parse_fits_wcs(filter(pair -> !startswith(pair.first, "A_") && !startswith(pair.first, "B_"), scamp_header))
    @test pixel_to_world(scamp, 500.0, 500.0) ≈ pixel_to_world(scamp_without_sip, 500.0, 500.0) atol=1e-13
    @test scamp_header == original_scamp_header

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
        code in ("COP", "COD", "COE", "COO", "BON") && (header["PV2_1"] = 45.0)
        world = pixel_to_world(parse_fits_wcs(header), 11.0, 22.0)
        tolerance = code == "CSC" ? 1e-8 : 5e-12
        @test world ≈ SVector(expected) atol=tolerance
    end
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
