# Audits inference and allocations of the scalar and matrix evaluation API over representative WCS pipelines:
# every FITS projection, distortions, spectral, time, tabular, and mixed-axis headers, the ASDF fixtures, and
# programmatic pipelines. Pass Roman or other ASDF files with a `roman.meta.wcs` or `wcs` node to include them.
#
#     julia --project=validation validation/audit_inference.jl [FILE.asdf...]
#
# A check fails when a call's return type is not inferred, when a warm scalar call allocates or throws, or when a
# matrix call allocates per column. Inverse checks use `quiet=true`, so convergence does not mask them, and are
# "n/a" for a WCS that has no analytical inverse and different pixel and world dimensions. The script exits with
# status 1 if any check fails.

using ASDF
using GWCS
using GWCS: AffineTransform, CallableLabelMapper, CallableTransform, CartesianToSpherical, CelestialFrame, CoordinateFrame, FITSTabularTable,
            Frame2D, NoInverse, Polynomial1D, Polynomial2D, RegionSelector, RemapAxes, Rotation3D, Scale, Shift, SphericalToCartesian, Step,
            WithInverseLabelMapper, compose, fix_inputs, input_dimension, output_dimension, parallel, parse_transform
using Logging
using StaticArrays
using Test

const FIXTURES = joinpath(@__DIR__, "..", "test", "fixtures")

@noinline allocated(f::F, args::Vararg{Any,N}) where {F,N} = @allocated f(args...)
@noinline allocations(f::F, args::Vararg{Any,N}) where {F,N} = @allocations f(args...)

function check(f, args...)
    try
        f(args...)
    catch error
        return "throws $(nameof(typeof(error)))"
    end
    inferred = try
        @inferred f(args...)
        true
    catch error
        error isa ErrorException || rethrow()
        false
    end
    bytes = allocated(f, args...)
    return inferred && bytes == 0 ? "ok" : !inferred && bytes == 0 ? "Any" : inferred ? "$(bytes)B" : "Any $(bytes)B"
end

# Allocations that grow with the number of columns mean a per-coordinate allocation in the batch path.
function per_column(f, column)
    small, large = repeat(column, 1, 10), repeat(column, 1, 1000)
    try
        f(small); f(large)
    catch error
        return "throws $(nameof(typeof(error)))"
    end
    # A per-column allocation would add at least 990; a few extra come from allocating a larger output.
    extra = allocations(f, large) - allocations(f, small)
    return extra < 10 ? "ok" : "$(round(extra / 990; digits=2))/col"
end

fits(header) = parse_fits_wcs(header; legacy=:quiet)
sky(code; extra...) = merge(Dict{String,Any}("NAXIS" => 2, "CTYPE1" => "RA---$code", "CTYPE2" => "DEC--$code", "CRPIX1" => 10.0, "CRPIX2" => 20.0,
    "CRVAL1" => 30.0, "CRVAL2" => 10.0, "CDELT1" => -0.1, "CDELT2" => 0.1), Dict{String,Any}(string(key) => value for (key, value) in extra))
load_asdf(path) = with_logger(() -> ASDF.load(path; extensions=true, validate_checksum=false), NullLogger())
wrap(transform) = GWCS.WCS(Step(CoordinateFrame{input_dimension(transform)}(name="input"), transform), Step(CoordinateFrame{output_dimension(transform)}(name="output")))

function cases(files)
    list = Pair{String,Any}[]
    projection_parameters = Dict("CEA" => (PV2_1=0.75,), "ZPN" => (PV2_1=1.0,), "COP" => (PV2_1=45.0,), "COD" => (PV2_1=45.0,), "COE" => (PV2_1=45.0,), "COO" => (PV2_1=45.0,), "BON" => (PV2_1=45.0,),
        "AZP" => (PV2_1=1.5, PV2_2=15.0), "SZP" => (PV2_1=1.5, PV2_2=15.0, PV2_3=80.0), "SIN" => (PV2_1=0.1, PV2_2=0.2), "AIR" => (PV2_1=75.8,), "CYP" => (PV2_1=1.0, PV2_2=1.0))
    for code in ("TAN", "SIN", "STG", "ARC", "ZEA", "AZP", "SZP", "AIR", "ZPN", "CAR", "CEA", "CYP", "MER", "SFL", "PAR", "MOL", "PCO", "AIT",
                 "COP", "COD", "COE", "COO", "BON", "TSC", "CSC", "QSC", "HPX", "XPH")
        push!(list, "FITS $code" => (fits(sky(code; get(projection_parameters, code, (;))...)), SVector(11.0, 22.0)))
    end
    sip = (A_ORDER=2, A_2_0=5e-6, A_1_1=-2e-6, A_0_2=3e-6, B_ORDER=2, B_2_0=-1e-6, B_1_1=4e-6, B_0_2=2e-6)
    reverse_sip = (AP_ORDER=2, AP_2_0=-5e-6, AP_1_1=2e-6, AP_0_2=-3e-6, BP_ORDER=2, BP_2_0=1e-6, BP_1_1=-4e-6, BP_0_2=-2e-6)
    push!(list, "FITS TAN-SIP" => (fits(sky("TAN-SIP"; sip..., reverse_sip...)), SVector(11.0, 22.0)))
    push!(list, "FITS TAN-SIP numerical" => (fits(sky("TAN-SIP"; sip...)), SVector(11.0, 22.0)))
    push!(list, "FITS TPV" => (fits(sky("TPV"; PV1_1=1.0, PV1_4=2e-3, PV1_5=-1e-3, PV1_7=5e-4, PV2_1=1.0, PV2_4=-1e-3, PV2_6=2e-3, PV2_8=-4e-4)), SVector(11.0, 22.0)))
    tpd = Dict{String,Any}("CPDIS1" => "TPD", "CPDIS2" => "TPD", "DP1.NAXES" => 2, "DP1.AXIS.1" => 1, "DP1.AXIS.2" => 2, "DP2.NAXES" => 2, "DP2.AXIS.1" => 2, "DP2.AXIS.2" => 1,
        "DP1.TPD.FWD.1" => 1.0, "DP1.TPD.FWD.4" => 1e-4, "DP1.TPD.FWD.5" => -2e-4, "DP2.TPD.FWD.1" => 1.0, "DP2.TPD.FWD.4" => -1e-4, "DP2.TPD.FWD.6" => 3e-4)
    push!(list, "FITS TAN + TPD distortion" => (fits(merge(sky("TAN"), tpd)), SVector(11.0, 22.0)))
    push!(list, "FITS ZPN + SIP-like TPD reverse" => (fits(merge(sky("ZPN"; PV2_1=1.0), merge(tpd, Dict("DP1.TPD.REV.1" => 1.0, "DP2.TPD.REV.1" => 1.0)))), SVector(11.0, 22.0)))
    linear = Dict{String,Any}("NAXIS" => 1, "CRPIX1" => 1.0, "CDELT1" => 1.0)
    push!(list, "FITS LINEAR" => (fits(merge(linear, Dict("CTYPE1" => "LINEAR", "CRVAL1" => 5.0))), SVector(2.0)))
    push!(list, "FITS FREQ-LOG" => (fits(merge(linear, Dict("CTYPE1" => "FREQ-LOG", "CUNIT1" => "Hz", "CRVAL1" => 1.4e9, "CDELT1" => 1e6))), SVector(2.0)))
    push!(list, "FITS WAVE-GRI" => (fits(merge(linear, Dict("CTYPE1" => "WAVE-GRI", "CUNIT1" => "m", "CRVAL1" => 5e-7, "CDELT1" => 1e-9, "PV1_0" => 3.16e5, "PV1_1" => 1.0, "PV1_2" => 13.9))), SVector(2.0)))
    push!(list, "FITS VELO-F2V" => (fits(merge(linear, Dict("CTYPE1" => "VELO-F2V", "CUNIT1" => "m/s", "CRVAL1" => 1000.0, "CDELT1" => 100.0, "RESTFRQ" => 1.420405751e9))), SVector(2.0)))
    push!(list, "FITS AWAV-V2A" => (fits(merge(linear, Dict("CTYPE1" => "AWAV-V2A", "CUNIT1" => "m", "CRVAL1" => 5e-7, "CDELT1" => 1e-10, "RESTFRQ" => 6e14))), SVector(2.0)))
    push!(list, "FITS TIME" => (fits(Dict("NAXIS" => 1, "CTYPE1" => "TIME", "TIMEUNIT" => "d", "TIMESYS" => "UTC", "MJDREFI" => 59000, "MJDREFF" => 0.25, "CRPIX1" => 1.0, "CRVAL1" => 1.0, "CDELT1" => 0.5)), SVector(2.0)))
    cube = Dict{String,Any}("NAXIS" => 4, "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CTYPE3" => "FREQ", "CTYPE4" => "STOKES", "CUNIT3" => "Hz",
        "CRPIX1" => 10.0, "CRPIX2" => 20.0, "CRPIX3" => 1.0, "CRPIX4" => 1.0, "CRVAL1" => 30.0, "CRVAL2" => 10.0, "CRVAL3" => 1.4e9, "CRVAL4" => 1.0,
        "CDELT1" => -0.1, "CDELT2" => 0.1, "CDELT3" => 1e6, "CDELT4" => 1.0)
    push!(list, "FITS RA/DEC/FREQ/STOKES cube" => (fits(cube), SVector(11.0, 22.0, 3.0, 1.0)))
    sip_cube = merge(cube, Dict{String,Any}("CTYPE1" => "RA---TAN-SIP", "CTYPE2" => "DEC--TAN-SIP"), Dict{String,Any}(string(key) => value for (key, value) in pairs(sip)))
    push!(list, "FITS TAN-SIP numerical cube" => (fits(sip_cube), SVector(11.0, 22.0, 3.0, 1.0)))
    interleaved = merge(cube, Dict("NAXIS" => 3, "CTYPE1" => "FREQ", "CUNIT1" => "Hz", "CTYPE3" => "RA---TAN", "CUNIT3" => "deg", "CRPIX1" => 1.0, "CRPIX3" => 10.0, "CRVAL1" => 1.4e9, "CRVAL3" => 30.0, "CDELT1" => 1e6, "CDELT3" => -0.1))
    foreach(key -> delete!(interleaved, key), ("CTYPE4", "CRPIX4", "CRVAL4", "CDELT4"))
    push!(list, "FITS FREQ/DEC/RA interleaved" => (fits(interleaved), SVector(2.0, 22.0, 11.0)))
    # Table probes sit inside the tables, at zero-based pixels.
    tab = Dict("NAXIS" => 1, "CTYPE1" => "FREQ-TAB", "CUNIT1" => "GHz", "CRPIX1" => 1.0, "CRVAL1" => 1.0, "CDELT1" => 1.0, "PS1_0" => "WCS-TABLE", "PS1_1" => "COORD")
    push!(list, "FITS FREQ-TAB" => (GWCS._parse_fits_wcs(GWCS._fits_header(tab); tabular_data=(FITSTabularTable((1,), reshape([10.0, 20.0, 40.0], 1, 3), ([1.0, 2.0, 3.0],)),)), SVector(1.5)))
    coupled = Dict("NAXIS" => 2, "CTYPE1" => "RA---TAB", "CTYPE2" => "DEC--TAB", "CRPIX1" => 1.0, "CRPIX2" => 1.0, "CRVAL1" => 1.0, "CRVAL2" => 1.0, "CDELT1" => 1.0, "CDELT2" => 1.0,
        "PS1_0" => "WCS-TABLE", "PS2_0" => "WCS-TABLE", "PS1_1" => "COORD", "PS2_1" => "COORD", "PV1_3" => 1, "PV2_3" => 2)
    coordinates = [axis == 1 ? 10i + j : i + 10j for axis in 1:2, i in 1:2, j in 1:2]
    push!(list, "FITS RA/DEC-TAB" => (GWCS._parse_fits_wcs(GWCS._fits_header(coupled); tabular_data=(FITSTabularTable((1, 2), Float64.(coordinates), ([1.0, 2.0], [1.0, 2.0])),)), SVector(0.5, 0.5)))

    push!(list, "ASDF roman_wcs.yaml" => (parse_wcs(read(joinpath(FIXTURES, "roman_wcs.yaml"), String)), SVector(11.0, 22.0)))
    push!(list, "ASDF asdf_parser_wcs.yaml" => (parse_wcs(read(joinpath(FIXTURES, "asdf_parser_wcs.yaml"), String)), SVector(11.0, 22.0)))
    mixed = parse_wcs(read(joinpath(FIXTURES, "mixed_frame_wcs.yaml"), String))
    push!(list, "ASDF mixed_frame_wcs.yaml" => (mixed, SVector(ntuple(i -> 1.0 + i, input_dimension(mixed.forward)))))
    quantities = with_logger(NullLogger()) do
        ASDF.load(joinpath(FIXTURES, "quantities.asdf"); extensions=true)
    end
    for name in ("mixed", "fiducial")
        push!(list, "ASDF quantities.asdf $name" => (parse_wcs(quantities["wcs"][name]), SVector(1.0, 2.0)))
    end
    projections = with_logger(NullLogger()) do
        ASDF.load(joinpath(FIXTURES, "projections.asdf"); extensions=true)
    end
    for name in ("car", "bon")
        push!(list, "ASDF projections.asdf $name" => (parse_wcs(projections["wcs"][name]), SVector(1.0, 2.0)))
    end
    boxes = GWCS._load_tagged_yaml(read(joinpath(FIXTURES, "asdf_bounding_boxes.yaml"), String))
    push!(list, "ASDF bounding box (ignored axis)" => (parse_wcs(boxes["simple_ignored"]), SVector(2.0, 3.0, 1.0)))
    push!(list, "ASDF compound bounding box" => (parse_wcs(boxes["compound"]), SVector(2.0, 3.0, 1.0)))
    selectors = GWCS._load_tagged_yaml(read(joinpath(FIXTURES, "asdf_selectors.yaml"), String))
    for (key, pixel) in (("array", SVector(0.0, 1.0)), ("dictionary", SVector(3.0, 8.0)), ("range", SVector(1.0, 8.0)), ("transform", SVector(99.0, 8.0)))
        push!(list, "ASDF $key region selector" => (wrap(parse_transform(selectors[key])), pixel))
    end
    arithmetic = GWCS._load_tagged_yaml(read(joinpath(FIXTURES, "asdf_arithmetic.yaml"), String))["arithmetic"]
    for key in sort(collect(keys(arithmetic)))
        push!(list, "ASDF arithmetic $key" => (wrap(parse_transform(arithmetic[key])), SVector(0.5)))
    end
    sky_correction = load_asdf(joinpath(FIXTURES, "sky_correction.asdf"))
    push!(list, "ASDF sky_correction.asdf (fallback)" => (parse_wcs(sky_correction["wcs"]), SVector(512.0, 512.0)))

    pixels, world = Frame2D(name="pixels"), CoordinateFrame{2}(name="world")
    rotation = parallel(Scale(1e-4), Scale(1e-4)) | SphericalToCartesian() | Rotation3D([0.01, -0.02, 30.0, -5.4, 83.8], "zyxyz") | CartesianToSpherical()
    push!(list, "sky rotation 2-3-3-2" => (GWCS.WCS(Step(pixels, rotation), Step(CelestialFrame(name="sky"))), SVector(400.0, 300.0)))
    pair = RemapAxes{2}((1, 2, 1, 2)) | parallel(Polynomial2D([0.1 0.2 0.3; 1.0 0.01 0.0; 0.02 0.0 0.0]), Polynomial2D([0.5 1.1 0.03; -0.2 0.0 0.0; 0.04 0.0 0.0]))
    push!(list, "polynomial pair distortion" => (GWCS.WCS(Step(pixels, pair | Shift(1.0) & Shift(2.0)), Step(world)), SVector(0.3, 0.4)))
    callable = CallableTransform{2,2}((x, y) -> (2x + y + 1, -x + 3y - 2); inverse=(u, v) -> ((3u - v - 5) / 7, (u + 2v + 3) / 7))
    push!(list, "callable (ForwardDiff Jacobian)" => (GWCS.WCS(Step(pixels, callable), Step(world)), SVector(3.0, 4.0)))
    push!(list, "elementwise arithmetic" => (wrap(Shift(1.0) * Scale(2.0)), SVector(3.0)))
    paired = WithInverseLabelMapper(CallableLabelMapper{2,Int}((x, y) -> x < 5 ? 1 : 2), CallableLabelMapper{2,Int}((x, y) -> x < 6 ? 1 : 2))
    regions = RegionSelector(paired, (1 => (Shift(1.0) & Shift(1.0)), 2 => (Shift(2.0) & Shift(2.0))))
    distortion = parallel(Polynomial1D([0.0, 1.0, 1e-3]), Polynomial1D([0.0, 1.0, -5e-4]))
    push!(list, "invertible selector + distortion (numerical)" => (GWCS.WCS(Step(pixels, regions | distortion), Step(world)), SVector(2.5, 3.5)))
    fixed = fix_inputs(GWCS.WCS(Step(CoordinateFrame{3}(name="cube", axes_names=("x", "y", "slice")), AffineTransform(SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0))), Step(CoordinateFrame{3}(name="out"))), "slice" => 2.0)
    push!(list, "fix_inputs" => (fixed, SVector(1.0, 2.0)))
    permuted = GWCS.WCS(Step(CoordinateFrame{3}(name="in", axes_order=(2, 3, 1)), AffineTransform(SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0))), Step(CoordinateFrame{3}(name="out", axes_order=(3, 1, 2))))
    push!(list, "permuted axes_order" => (permuted, SVector(1.0, 2.0, 3.0)))

    for file in files
        tree = load_asdf(file)
        node = haskey(tree.metadata, "roman") ? tree["roman"]["meta"]["wcs"] : tree["wcs"]
        wcs = parse_wcs(node)
        center = wcs.pixel_shape === nothing ? SVector(ntuple(_ -> 10.0, input_dimension(wcs.forward))) : SVector(map(n -> (n - 1) / 2, wcs.pixel_shape))
        push!(list, basename(file) => (wcs, center))
        strip(transform) = transform isa GWCS.WithInverse ? transform.forward : transform
        if wcs.forward isa GWCS.ComposedTransform
            numerical = GWCS.WCS(Step(input_frame(wcs), compose(map(strip, wcs.forward.transforms)...)), Step(output_frame(wcs)); bounding_box=wcs.bounding_box)
            numerical.backward isa NoInverse && push!(list, basename(file) * " (numerical inverse)" => (numerical, center))
        end
    end
    return list
end

function audit(name, wcs, pixel)
    world = pixel_to_world(wcs, pixel; with_bounding_box=false)
    all(isfinite, world) || return (name, fill("no finite world", 10))
    # Without an analytical inverse, a WCS whose pixel and world dimensions differ cannot be inverted at all.
    invertible = !(wcs.backward isa NoInverse) || length(pixel) == length(world)
    inverse(f) = invertible ? f() : "n/a"
    checks = String[
        check(p -> pixel_to_world(wcs, p), pixel),
        check(p -> pixel_to_world(wcs, p...), Tuple(pixel)),
        inverse(() -> check(w -> world_to_pixel(wcs, w; quiet=true), world)),
        check(p -> pixel_to_world_jacobian(wcs, p; with_bounding_box=false), pixel),
        inverse(() -> check(w -> world_to_pixel_jacobian(wcs, w; with_bounding_box=false, quiet=true), world)),
        check(p -> in_domain(wcs, p), pixel),
        inverse(() -> check(w -> in_image(wcs, w), world)),
        check(p -> pixel_to_world(wcs, p), pixel .+ 1e7),
        per_column(m -> pixel_to_world(wcs, m), Vector(pixel)),
        inverse(() -> per_column(m -> world_to_pixel(wcs, m; quiet=true), Vector(world))),
    ]
    return (name, checks)
end

const COLUMNS = ("pixel_to_world", "p2w(x, y...)", "world_to_pixel", "p2w jacobian", "w2p jacobian", "in_domain", "in_image", "p2w far/fill", "p2w matrix", "w2p matrix")

function main(files)
    rows = [audit(name, case...) for (name, case) in cases(files)]
    width = maximum(row -> length(row[1]), rows) + 2
    println(rpad("case", width), join(rpad.(COLUMNS, 18)))
    failures = 0
    for (name, checks) in rows
        failures += count(check -> check ∉ ("ok", "n/a"), checks)
        println(rpad(name, width), join(rpad.(checks, 18)))
    end
    println("\n", length(rows), " cases, ", failures, " failing checks (\"Any\": return type not inferred; \"NB\": bytes allocated per warm call; \"/col\": allocations per matrix column)")
    failures == 0 || exit(1)
end

main(ARGS)
