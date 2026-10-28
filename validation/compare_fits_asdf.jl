using ASDF
using GWCS
using JSON

include("python.jl")

# GWCS.jl writes FITS-derived WCSs, with their frames and polynomial distortions (SIP, TPV, and Paper IV TPD) in
# standard ASDF form; Python GWCS must read the same celestial frames and reference times that astropy derives from the
# headers, and the same coordinates, which astropy evaluates with the distortions.
const HEADERS = Dict(
    "icrs" => Dict("CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CRVAL1" => 30.0, "CRVAL2" => 40.0),
    "fk5" => Dict("CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CRVAL1" => 30.0, "CRVAL2" => 40.0, "RADESYS" => "FK5", "EQUINOX" => 1984.5),
    "fk4" => Dict("CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CRVAL1" => 30.0, "CRVAL2" => 40.0, "EQUINOX" => 1950.0),
    "fk4_no_e" => Dict("CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CRVAL1" => 30.0, "CRVAL2" => 40.0, "RADESYS" => "FK4-NO-E", "EQUINOX" => 1975.0),
    "galactic" => Dict("CTYPE1" => "GLON-CAR", "CTYPE2" => "GLAT-CAR", "CRVAL1" => 120.0, "CRVAL2" => -10.0),
    "dec_first" => Dict("CTYPE1" => "DEC--TAN", "CTYPE2" => "RA---TAN", "CRVAL1" => 40.0, "CRVAL2" => 30.0, "RADESYS" => "FK5", "EQUINOX" => 2000.0),
    "cube" => Dict("CTYPE1" => "RA---TAN", "CTYPE2" => "WAVE", "CTYPE3" => "DEC--TAN", "CRVAL1" => 30.0, "CRVAL2" => 5e-7, "CRVAL3" => 40.0, "CUNIT2" => "m"),
    "time_topocentric" => Dict("CTYPE1" => "TIME", "CRVAL1" => 10.0, "CTYPE2" => "X", "CRVAL2" => 1.0, "MJDREFI" => 59000, "MJDREFF" => 0.25, "TIMESYS" => "TT", "TREFPOS" => "TOPOCENTER",
        "OBSGEO-X" => 5.1e6, "OBSGEO-Y" => -1.2e6, "OBSGEO-Z" => 3.6e6),
    "time_alone" => Dict("CTYPE1" => "TIME", "CRVAL1" => 10.0, "MJDREF" => 59000.0, "TIMESYS" => "UTC"),
    "time_gps" => Dict("CTYPE1" => "TIME", "CRVAL1" => 10.0, "CTYPE2" => "X", "CRVAL2" => 1.0, "MJDREF" => 59000.0, "TIMESYS" => "GPS"),
    "sip" => Dict("CTYPE1" => "RA---TAN-SIP", "CTYPE2" => "DEC--TAN-SIP", "CRVAL1" => 83.8, "CRVAL2" => -5.4, "CRPIX1" => 512.0, "CRPIX2" => 512.0,
        "A_ORDER" => 3, "A_2_0" => 2e-6, "A_1_1" => -1e-6, "A_0_3" => 3e-9, "B_ORDER" => 3, "B_0_2" => -2e-6, "B_2_1" => 1e-9,
        "AP_ORDER" => 2, "AP_2_0" => -2e-6, "BP_ORDER" => 2, "BP_0_2" => 2e-6),
    "tpv" => Dict("CTYPE1" => "RA---TPV", "CTYPE2" => "DEC--TPV", "CRVAL1" => 83.8, "CRVAL2" => -5.4, "CRPIX1" => 512.0, "CRPIX2" => 512.0,
        "PV1_1" => 1.0, "PV1_4" => 0.01, "PV1_5" => -0.003, "PV2_1" => 1.0, "PV2_4" => -0.02, "PV2_6" => 0.004),
    "prior_tpd" => Dict("CTYPE1" => "LINEAR", "CTYPE2" => "LINEAR", "CRVAL1" => 5.0, "CRVAL2" => -7.0, "CRPIX1" => 10.0, "CRPIX2" => 20.0,
        "CPDIS1" => "TPD", "CPDIS2" => "TPD", "DP1.NAXES" => 2, "DP1.AXIS.1" => 1, "DP1.AXIS.2" => 2, "DP1.OFFSET.1" => 1.0, "DP1.OFFSET.2" => -2.0,
        "DP1.SCALE.1" => 0.5, "DP1.SCALE.2" => 2.0, "DP1.TPD.FWD.0" => 0.01, "DP1.TPD.FWD.1" => 1.02, "DP1.TPD.FWD.2" => -0.03, "DP1.TPD.FWD.4" => 0.0007,
        "DP2.NAXES" => 2, "DP2.AXIS.1" => 2, "DP2.AXIS.2" => 1, "DP2.TPD.FWD.4" => -0.002, "DP2.TPD.FWD.5" => 0.0003),
    "sequent_tpd" => Dict("CTYPE1" => "LINEAR", "CTYPE2" => "LINEAR", "CRVAL1" => 5.0, "CRVAL2" => -7.0, "CRPIX1" => 10.0, "CRPIX2" => 20.0,
        "CQDIS1" => "TPD", "CQDIS2" => "TPD", "DQ1.NAXES" => 2, "DQ1.AXIS.1" => 1, "DQ1.AXIS.2" => 2, "DQ1.OFFSET.1" => 1.0, "DQ1.OFFSET.2" => -2.0,
        "DQ1.SCALE.1" => 0.5, "DQ1.SCALE.2" => 2.0, "DQ1.DOCORR" => 0, "DQ1.AUX.1.COEFF.0" => 0.1, "DQ1.AUX.1.COEFF.1" => 1.1, "DQ1.AUX.1.COEFF.2" => -0.2,
        "DQ1.AUX.2.COEFF.0" => -0.3, "DQ1.AUX.2.COEFF.1" => 0.4, "DQ1.AUX.2.COEFF.2" => 0.9, "DQ1.TPD.FWD.0" => 0.01, "DQ1.TPD.FWD.1" => 1.02,
        "DQ1.TPD.FWD.4" => 0.0007, "DQ2.NAXES" => 2, "DQ2.AXIS.1" => 2, "DQ2.AXIS.2" => 1, "DQ2.TPD.FWD.4" => -0.002),
    "time_dateref" => Dict("CTYPE1" => "TDB", "CRVAL1" => 10.0, "CTYPE2" => "X", "CRVAL2" => 1.0, "DATEREF" => "2020-05-31T06:00:00", "TREFPOS" => "GEOCENTER"),
)

function header(cards)
    count = maximum(parse(Int, match(r"\d+$", key).match) for key in keys(cards) if startswith(key, "CTYPE"))
    return merge(Dict{String,Any}("NAXIS" => count, ("CRPIX$i" => 1.0 for i in 1:count)..., ("CDELT$i" => i == 2 && get(cards, "CTYPE2", "") == "WAVE" ? 1e-9 : 1e-3 for i in 1:count)...), cards)
end

function main()
    headers = Dict(name => header(cards) for (name, cards) in HEADERS)
    wcss = Dict(name => parse_fits_wcs(cards) for (name, cards) in headers)
    mktempdir() do directory
        filename, header_file = joinpath(directory, "fits_asdf.asdf"), joinpath(directory, "headers.json")
        write_asdf(filename, Dict("wcs" => wcss))
        write(header_file, JSON.json(headers))
        records = JSON.parse(read(`$(validation_python()) $(joinpath(@__DIR__, "python_fits_asdf_reference.py")) $filename $header_file`, String))
        for record in sort(records; by=record -> record["name"])
            name = record["name"]
            # GWCS.jl's own coordinates, which the WCSLIB comparisons check, for every header; astropy's where it reads it.
            julia_error = maximum(abs, pixel_to_world(wcss[name], Float64.(record["pixel"])) .- Float64.(record["values"]))
            astropy = haskey(record, "astropy_error") ? "astropy rejects: $(record["astropy_error"])" : "$(rpad(record["frame"], 48)) astropy error $(record["astropy_error_value"])"
            println(rpad(name, 18), "GWCS.jl error ", rpad(julia_error, 24), astropy)
            julia_error <= 1e-9 || error("$name: Python reads coordinates that differ from GWCS.jl's by $julia_error")
            haskey(record, "astropy_error") && continue
            record["frame_matches"] || error("$name: Python reads $(record["frame"]), astropy gives $(record["astropy_frame"])")
            record["astropy_error_value"] <= 1e-9 || error("$name: coordinates differ from astropy by $(record["astropy_error_value"])")
        end
    end
end

main()
