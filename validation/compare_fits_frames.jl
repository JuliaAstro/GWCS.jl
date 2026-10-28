using ASDF
using GWCS
using JSON

include("python.jl")

# GWCS.jl writes FITS-derived WCSs whose frames have an ASDF form; Python GWCS must read the same celestial frames and
# reference times that astropy derives from the headers, and the same coordinates.
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
        filename, header_file = joinpath(directory, "fits_frames.asdf"), joinpath(directory, "headers.json")
        write_asdf(filename, Dict("wcs" => wcss))
        write(header_file, JSON.json(headers))
        records = JSON.parse(read(`$(validation_python()) $(joinpath(@__DIR__, "python_fits_frames_reference.py")) $filename $header_file`, String))
        for record in sort(records; by=record -> record["name"])
            println(rpad(record["name"], 18), rpad(record["frame"], 48), " value error ", record["value_error"])
            record["frame_matches"] || error("$(record["name"]): Python reads $(record["frame"]), astropy gives $(record["astropy_frame"])")
            record["value_error"] <= 1e-9 || error("$(record["name"]): coordinates differ from astropy by $(record["value_error"])")
        end
    end
end

main()
