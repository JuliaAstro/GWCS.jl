"""
First-call latency of constructing and evaluating a metadata-defined WCS, which is dominated by compilation and so
invisible to BenchmarkTools. Run it in a fresh process:

    julia --project=benchmark benchmark/latency.jl [FILE.yaml]

The default is the Roman fixture in `test/fixtures/roman_wcs.yaml`.
"""

stage(label, f) = (time = @elapsed(result = Base.invokelatest(f)); println(rpad(label, 24), round(time; digits=3), " s"); result)

stage("using GWCS", () -> @eval(using GWCS, StaticArrays))
path = isempty(ARGS) ? joinpath(@__DIR__, "..", "test", "fixtures", "roman_wcs.yaml") : only(ARGS)
tree = stage("load metadata", () -> GWCS._load_tagged_yaml(read(path, String)))
wcs = stage("first parse_wcs", () -> parse_wcs(tree))
pixel = SVector(ntuple(_ -> 1.0, GWCS.input_dimension(wcs.forward)))
world = stage("first pixel_to_world", () -> pixel_to_world(wcs, pixel))
stage("first world_to_pixel", () -> world_to_pixel(wcs, world))
stage("second parse_wcs", () -> parse_wcs(tree))
