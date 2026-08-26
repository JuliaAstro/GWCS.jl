using ASDF
using GWCS

length(ARGS) == 1 || error("usage: julia --project=validation validation/load_roman.jl FILE")

asdf = ASDF.load(only(ARGS); extensions=true, validate_checksum=false)
wcs = parse_wcs(asdf["roman"]["meta"]["wcs"])
println("frames: ", join(available_frames(wcs), " -> "))
println("pixel (2043.5, 2043.5): ", pixel_to_world(wcs, 2043.5, 2043.5))
