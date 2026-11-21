using ASDF
using GWCS
using GWCS: celestial_system, with_celestial_system
using JSON
using StaticArrays

include("python.jl")

# WCSs converted with `with_celestial_system` and written to ASDF must give the same coordinates when Python GWCS reads
# them, in the astropy frame of the new system, and astropy must give those coordinates when it evaluates the original
# (a FITS header, or the original WCS as Python GWCS reads it) and converts them with `SkyCoord.transform_to`. Replaced
# output frames are fused into the original's final rotation where it has one; kept frames are followed by their own
# rotation step.
function main()
    tan(ctype1, ctype2, extra...) = Dict{String,Any}("CTYPE1" => ctype1, "CTYPE2" => ctype2, "CRPIX1" => 50.0, "CRPIX2" => 40.0, "CRVAL1" => 120.0,
        "CRVAL2" => 20.0, "CDELT1" => -1e-3, "CDELT2" => 1e-3, "PC1_2" => 0.1, extra...)
    roman = parse_wcs(read(joinpath(@__DIR__, "..", "test", "fixtures", "roman_wcs.yaml"), String))
    icrs, galactic, supergalactic = ("icrs", 0.0), ("galactic", 0.0), ("supergalactic", 0.0)
    # name => (original, its system, target system, target equinox, keep_frame)
    cases = Dict(
        "galactic_to_icrs" => (tan("GLON-TAN", "GLAT-TAN"), galactic, :icrs, nothing, false),
        "galactic_to_icrs_kept" => (tan("GLON-TAN", "GLAT-TAN"), galactic, :icrs, nothing, true),
        "icrs_to_galactic" => (tan("RA---TAN", "DEC--TAN"), icrs, :galactic, nothing, false),
        "fk5_1975_to_icrs" => (tan("RA---TAN", "DEC--TAN", "RADESYS" => "FK5", "EQUINOX" => 1975.0), ("fk5", 1975.0), :icrs, nothing, false),
        "icrs_to_fk5_2025" => (tan("RA---TAN", "DEC--TAN"), icrs, :fk5, 2025.5, false),
        "galactic_to_fk5_1975_kept" => (tan("GLON-TAN", "GLAT-TAN"), galactic, :fk5, 1975.0, true),
        "reversed_axes" => (tan("DEC--TAN", "RA---TAN", "CRVAL1" => 20.0, "CRVAL2" => 120.0), icrs, :galactic, nothing, false),
        "supergalactic_to_icrs" => (tan("SLON-TAN", "SLAT-TAN"), supergalactic, :icrs, nothing, false),
        "roman_to_galactic" => (roman, icrs, :galactic, nothing, false),
        "roman_to_galactic_kept" => (roman, icrs, :galactic, nothing, true),
    )
    pixels = [[0.0, 0.0], [10.0, 20.0], [19.0, 29.0]]
    convert(original) = original isa WCS ? original : parse_fits_wcs(original)
    wcs = Dict(name => with_celestial_system(convert(original), target; keep_frame, (equinox === nothing ? (;) : (; equinox))...) for (name, (original, _, target, equinox, keep_frame)) in cases)
    systems(system, equinox) = Dict("system" => system, "equinox" => equinox)
    tree = Dict("wcs" => wcs, "pixels" => pixels,
        "headers" => Dict(name => original for (name, (original, _...)) in cases if !(original isa WCS)),
        "originals" => Dict(name => original for (name, (original, _...)) in cases if original isa WCS),
        "sources" => Dict(name => systems(source...) for (name, (_, source, _...)) in cases),
        "targets" => Dict(name => systems(String(target), something(equinox, 0.0)) for (name, (_, _, target, equinox, _)) in cases))
    max_gwcs, max_astropy, max_inverse = 0.0, 0.0, 0.0
    mktempdir() do directory
        filename = joinpath(directory, "converted.asdf")
        ASDF.write_file(filename, tree)
        records = JSON.parse(read(`$(validation_python()) $(joinpath(@__DIR__, "python_celestial_conversion_reference.py")) $filename`, String))
        length(records) == length(cases) || error("Python read $(length(records)) of $(length(cases)) WCSs")
        for record in sort(records; by=record -> record["wcs"])
            name = record["wcs"]
            system = celestial_system(GWCS.celestial_frame(wcs[name]))
            record["frame"] == String(system.system) || error("$name: Python reads a $(record["frame"]) frame, not $(system.system)")
            record["inverse"] === nothing && error("$name: Python GWCS has no analytical inverse")
            max_inverse = max(max_inverse, maximum(k -> maximum(abs, record["inverse"][k] .- pixels[k]), eachindex(pixels)))
            for (pixel, from_gwcs, from_astropy) in zip(pixels, record["gwcs"], record["astropy"])
                world = pixel_to_world(wcs[name], SVector{2}(pixel))
                max_gwcs = max(max_gwcs, maximum(abs, world .- from_gwcs))
                max_astropy = max(max_astropy, maximum(abs, world .- from_astropy))
            end
            println(rpad(name, 26), " Python reads ", record["frame"], ", ", length(wcs[name].steps), " frames")
        end
    end
    println("maximum difference from Python GWCS: ", max_gwcs, " deg; from astropy: ", max_astropy, " deg")
    println("maximum pixel error of Python GWCS's analytical inverse: ", max_inverse)
    max(max_gwcs, max_astropy) <= 1e-9 || error("celestial conversion comparison exceeded 1e-9 deg")
    max_inverse <= 1e-6 || error("Python GWCS's inverse exceeded 1e-6 pixels")
end

main()
