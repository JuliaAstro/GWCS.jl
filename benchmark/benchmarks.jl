"""
Comparative scalar benchmarks for equivalent FITSWCS.jl and GWCS.jl transforms.

Usage:
    julia --project=benchmark benchmark/benchmarks.jl
"""

using BenchmarkTools
using FITSWCS
using GWCS
using PrettyTables
using StaticArrays

function flatten_results(group)
    flat = Dict{String,Any}()
    _flatten_results!(flat, group, "")
    return flat
end

function _flatten_results!(flat, group, prefix)
    for (key, value) in group
        name = isempty(prefix) ? string(key) : "$prefix/$key"
        value isa BenchmarkGroup ? _flatten_results!(flat, value, name) : flat[name] = value
    end
    return flat
end

function show_benchmarks(results)
    sorted = sort(collect(flatten_results(results)); by=first)
    trials = Dict(sorted)
    ratios = map(sorted) do (name, value)
        endswith(name, "/GWCS") || return ""
        fits_name = replace(name, r"/GWCS$" => "/FITSWCS")
        haskey(trials, fits_name) || return ""
        return string(round(median(value).time / median(trials[fits_name]).time; digits=2), "x")
    end
    times = [BenchmarkTools.prettytime(median(value).time) for (key, value) in sorted]
    memory = [BenchmarkTools.prettymemory(median(value).memory) for (key, value) in sorted]
    allocations = [median(value).allocs for (key, value) in sorted]
    data = hcat(first.(sorted), times, memory, allocations, ratios)
    pretty_table(data; column_labels=["Benchmark", "Median Time", "Memory", "Allocs", "GWCS/FITSWCS"], alignment=[:l, :r, :r, :r, :r])
end

const PIXEL_FRAME = Frame2D(name="detector", axes_names=("x", "y"), unit=("pixel", "pixel"))
const WORLD_FRAME = CelestialFrame(name="world")
const LINEAR_FRAME = CoordinateFrame{2}(name="linear", axes_names=("x", "y"))

function affine_wcs(crpix, crval, matrix)
    transform = AffineTransform(matrix, crval - matrix * crpix)
    return GWCS.WCS(Step(PIXEL_FRAME, transform), Step(LINEAR_FRAME))
end

function tan_wcs(fits)
    linear = AffineTransform(fits.cd, -(fits.cd * fits.crpix))
    transform = compose(linear, TangentPlaneToNative(), NativeToCelestial(fits.alpha_p, fits.delta_p, fits.lonpole))
    return GWCS.WCS(Step(PIXEL_FRAME, transform), Step(WORLD_FRAME))
end

function sip_polynomial(crpix, coefficients)
    values = copy(coefficients)
    values[1, 1] += crpix
    return Polynomial2D(values)
end

function sip_distortion(crpix, a, b, ap, bp)
    center = parallel(Shift(-crpix[1]), Shift(-crpix[2]))
    duplicate = RemapAxes{2}((1, 2, 1, 2))
    forward = compose(center, duplicate, parallel(sip_polynomial(crpix[1], a + [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 0.0 0.0]), sip_polynomial(crpix[2], b + [0.0 1.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0])))
    backward = compose(center, duplicate, parallel(sip_polynomial(crpix[1], ap + [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 0.0 0.0]), sip_polynomial(crpix[2], bp + [0.0 1.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0])))
    return WithInverse(forward, backward)
end

function tan_sip_wcs(fits, a, b, ap, bp)
    distortion = sip_distortion(fits.crpix, a, b, ap, bp)
    linear = AffineTransform(fits.cd, -(fits.cd * fits.crpix))
    transform = compose(distortion, linear, TangentPlaneToNative(), NativeToCelestial(fits.alpha_p, fits.delta_p, fits.lonpole))
    return GWCS.WCS(Step(PIXEL_FRAME, transform), Step(WORLD_FRAME))
end

const LINEAR_HEADER = Dict(
    "NAXIS" => 2,
    "CTYPE1" => "LINEAR", "CTYPE2" => "LINEAR",
    "CRPIX1" => 512.0, "CRPIX2" => 512.0,
    "CRVAL1" => 12.0, "CRVAL2" => -8.0,
    "CD1_1" => 0.02, "CD1_2" => 0.001,
    "CD2_1" => -0.002, "CD2_2" => 0.03
)

const TAN_HEADER = Dict(
    "NAXIS" => 2,
    "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN",
    "CRPIX1" => 512.0, "CRPIX2" => 512.0,
    "CRVAL1" => 83.8221, "CRVAL2" => -5.3911,
    "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4
)

const SIP_HEADER = merge(copy(TAN_HEADER), Dict(
    "CTYPE1" => "RA---TAN-SIP", "CTYPE2" => "DEC--TAN-SIP",
    "A_ORDER" => 2,
    "A_2_0" => 5e-7, "A_1_1" => -2e-7, "A_0_2" => 3e-7,
    "B_ORDER" => 2,
    "B_2_0" => -1e-7, "B_1_1" => 4e-7, "B_0_2" => 2e-7,
    "AP_ORDER" => 2,
    "AP_2_0" => -5e-7, "AP_1_1" => 2e-7, "AP_0_2" => -3e-7,
    "BP_ORDER" => 2,
    "BP_2_0" => 1e-7, "BP_1_1" => -4e-7, "BP_0_2" => -2e-7
))

const FITS_LINEAR = FITSWCS.WCS(LINEAR_HEADER)
const FITS_TAN = FITSWCS.WCS(TAN_HEADER)
const FITS_SIP = FITSWCS.WCS(SIP_HEADER)

const GWCS_LINEAR = affine_wcs(FITS_LINEAR.crpix, FITS_LINEAR.crval, FITS_LINEAR.cd)
const GWCS_TAN = tan_wcs(FITS_TAN)
const SIP_A = FITS_SIP.pipeline.sip.a
const SIP_B = FITS_SIP.pipeline.sip.b
const SIP_AP = FITS_SIP.pipeline.sip.ap
const SIP_BP = FITS_SIP.pipeline.sip.bp
const GWCS_SIP = tan_sip_wcs(FITS_SIP, SIP_A, SIP_B, SIP_AP, SIP_BP)

callable_forward(x, y) = (2x + y + 1, -x + 3y - 2)
callable_inverse(u, v) = ((3u - v - 5) / 7, (u + 2v + 3) / 7)
callable_jacobian(x, y) = SMatrix{2,2}(2, -1, 1, 3)
arithmetic_forward(x) = (x + 1) * (2x)
arithmetic_jacobian(x) = SMatrix{1,1}(4x + 2)

const NATIVE_CALLABLE_REFERENCE = AffineTransform(SMatrix{2,2}(2.0, -1.0, 1.0, 3.0), SVector(1.0, -2.0))
const CALLABLE_REFERENCE = CallableTransform{2,2}(callable_forward; inverse=callable_inverse, jacobian=callable_jacobian)
const CALLABLE_AD_REFERENCE = CallableTransform{2,2}(callable_forward; inverse=callable_inverse)
const ARITHMETIC_REFERENCE = Shift(1.0) * Scale(2.0)
const ARITHMETIC_CALLABLE = CallableTransform{1,1}(arithmetic_forward; jacobian=arithmetic_jacobian)
const NATIVE_CALLABLE_INVERSE = inverse_transform(NATIVE_CALLABLE_REFERENCE)
const CALLABLE_INVERSE = inverse_transform(CALLABLE_REFERENCE)
const FIXED_REFERENCE = fix_inputs(NATIVE_CALLABLE_REFERENCE, 2 => 300.0)
const FIXED_INVERSE = inverse_transform(FIXED_REFERENCE)
const CONSTANT_AXIS = parallel(IdentityTransform{2}(), ConstantTransform(42.0))
const ORDERED_3D_FRAME = CoordinateFrame{3}(name="ordered")
const PERMUTED_INPUT_FRAME = CoordinateFrame{3}(name="permuted_input", axes_order=(2, 3, 1))
const PERMUTED_OUTPUT_FRAME = CoordinateFrame{3}(name="permuted_output", axes_order=(3, 1, 2))
const AXIS_ORDER_TRANSFORM = AffineTransform(SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0))
const ORDERED_3D_WCS = GWCS.WCS(Step(ORDERED_3D_FRAME, AXIS_ORDER_TRANSFORM), Step(CoordinateFrame{3}(name="ordered_output")))
const PERMUTED_3D_WCS = GWCS.WCS(Step(PERMUTED_INPUT_FRAME, AXIS_ORDER_TRANSFORM), Step(PERMUTED_OUTPUT_FRAME))
const DOMAIN_INPUT_FRAME = CoordinateFrame{3}(name="domain_input", axes_names=("x", "y", "slice"))
const DOMAIN_OUTPUT_FRAME = CoordinateFrame{3}(name="domain_output")
const SIMPLE_DOMAIN_BOX = BoundingBox{3}(((0.0, 5.0), (0.0, 5.0)), (1, 2))
const COMPOUND_DOMAIN_BOX = CompoundBoundingBox(Dict((1,) => SIMPLE_DOMAIN_BOX,
    (2,) => BoundingBox{3}(((10.0, 15.0), (10.0, 15.0)), (1, 2))), (3,))
const SIMPLE_DOMAIN_WCS = GWCS.WCS(Step(DOMAIN_INPUT_FRAME, IdentityTransform{3}()), Step(DOMAIN_OUTPUT_FRAME); bounding_box=SIMPLE_DOMAIN_BOX)
const COMPOUND_DOMAIN_WCS = GWCS.WCS(Step(DOMAIN_INPUT_FRAME, IdentityTransform{3}()), Step(DOMAIN_OUTPUT_FRAME); bounding_box=COMPOUND_DOMAIN_BOX)

const PIXEL = SVector(400.0, 300.0)
const FIXED_PIXEL = SVector(400.0)
const ARITHMETIC_PIXEL = SVector(3.0)
const PIXEL_3D = SVector(1.0, 2.0, 3.0)
const DOMAIN_INSIDE = SVector(2.0, 3.0, 1.0)
const DOMAIN_OUTSIDE = SVector(7.0, 3.0, 1.0)
const DOMAIN_BATCH = [2.0 7.0 12.0; 3.0 3.0 12.0; 1.0 1.0 2.0]
const CALLABLE_WORLD = CALLABLE_REFERENCE(PIXEL)
const FIXED_WORLD = FIXED_REFERENCE(FIXED_PIXEL)
const ORDERED_WORLD_3D = GWCS.pixel_to_world(ORDERED_3D_WCS, PIXEL_3D)
const PERMUTED_WORLD_3D = GWCS.pixel_to_world(PERMUTED_3D_WCS, PIXEL_3D)
const FITS_LINEAR_WORLD = FITSWCS.pixel_to_world(FITS_LINEAR, PIXEL)
const FITS_TAN_WORLD = FITSWCS.pixel_to_world(FITS_TAN, PIXEL)
const FITS_SIP_WORLD = FITSWCS.pixel_to_world(FITS_SIP, PIXEL)
const GWCS_LINEAR_WORLD = GWCS.pixel_to_world(GWCS_LINEAR, PIXEL)
const GWCS_TAN_WORLD = GWCS.pixel_to_world(GWCS_TAN, PIXEL)
const GWCS_SIP_WORLD = GWCS.pixel_to_world(GWCS_SIP, PIXEL)

@assert GWCS_LINEAR_WORLD ≈ FITS_LINEAR_WORLD rtol=1e-12
@assert GWCS_TAN_WORLD ≈ FITS_TAN_WORLD atol=1e-11
@assert GWCS_SIP_WORLD ≈ FITS_SIP_WORLD atol=1e-11
@assert GWCS.world_to_pixel(GWCS_LINEAR, FITS_LINEAR_WORLD) ≈ FITSWCS.world_to_pixel(FITS_LINEAR, FITS_LINEAR_WORLD) atol=1e-10
@assert GWCS.world_to_pixel(GWCS_TAN, FITS_TAN_WORLD) ≈ FITSWCS.world_to_pixel(FITS_TAN, FITS_TAN_WORLD) atol=1e-9
@assert GWCS.world_to_pixel(GWCS_SIP, FITS_SIP_WORLD) ≈ FITSWCS.world_to_pixel(FITS_SIP, FITS_SIP_WORLD) atol=1e-6

const SUITE = BenchmarkGroup()

for operation in ("pixel_to_world", "world_to_pixel", "construction", "callable", "arithmetic", "fixed_inputs", "axis_order", "domains")
    SUITE[operation] = BenchmarkGroup()
end

for (name, fits, gwcs, fits_world, gwcs_world) in (
    ("affine", FITS_LINEAR, GWCS_LINEAR, FITS_LINEAR_WORLD, GWCS_LINEAR_WORLD),
    ("TAN", FITS_TAN, GWCS_TAN, FITS_TAN_WORLD, GWCS_TAN_WORLD),
    ("TAN-SIP", FITS_SIP, GWCS_SIP, FITS_SIP_WORLD, GWCS_SIP_WORLD)
)
    SUITE["pixel_to_world"]["$name/FITSWCS"] = @benchmarkable FITSWCS.pixel_to_world($fits, $PIXEL) evals=100
    SUITE["pixel_to_world"]["$name/GWCS"] = @benchmarkable GWCS.pixel_to_world($gwcs, $PIXEL) evals=100
    SUITE["world_to_pixel"]["$name/FITSWCS"] = @benchmarkable FITSWCS.world_to_pixel($fits, $fits_world) evals=100
    SUITE["world_to_pixel"]["$name/GWCS"] = @benchmarkable GWCS.world_to_pixel($gwcs, $gwcs_world) evals=100
end

SUITE["construction"]["affine/GWCS"] = @benchmarkable affine_wcs($(FITS_LINEAR.crpix), $(FITS_LINEAR.crval), $(FITS_LINEAR.cd)) evals=1
SUITE["construction"]["TAN/GWCS"] = @benchmarkable tan_wcs($FITS_TAN) evals=1
SUITE["construction"]["TAN-SIP/GWCS"] = @benchmarkable tan_sip_wcs($FITS_SIP, $SIP_A, $SIP_B, $SIP_AP, $SIP_BP) evals=1

SUITE["callable"]["forward/native"] = @benchmarkable $NATIVE_CALLABLE_REFERENCE($PIXEL) evals=100
SUITE["callable"]["forward/wrapper"] = @benchmarkable $CALLABLE_REFERENCE($PIXEL) evals=100
SUITE["callable"]["forward/raw"] = @benchmarkable callable_forward($PIXEL[1], $PIXEL[2]) evals=100
SUITE["callable"]["inverse/native"] = @benchmarkable $NATIVE_CALLABLE_INVERSE($CALLABLE_WORLD) evals=100
SUITE["callable"]["inverse/wrapper"] = @benchmarkable $CALLABLE_INVERSE($CALLABLE_WORLD) evals=100
SUITE["callable"]["jacobian/analytical"] = @benchmarkable transform_jacobian($CALLABLE_REFERENCE, $PIXEL) evals=100
SUITE["callable"]["jacobian/ForwardDiff"] = @benchmarkable transform_jacobian($CALLABLE_AD_REFERENCE, $PIXEL) evals=100

SUITE["arithmetic"]["forward/structured"] = @benchmarkable $ARITHMETIC_REFERENCE($ARITHMETIC_PIXEL) evals=100
SUITE["arithmetic"]["forward/callable"] = @benchmarkable $ARITHMETIC_CALLABLE($ARITHMETIC_PIXEL) evals=100
SUITE["arithmetic"]["forward/raw"] = @benchmarkable arithmetic_forward($ARITHMETIC_PIXEL[1]) evals=100
SUITE["arithmetic"]["jacobian/structured"] = @benchmarkable transform_jacobian($ARITHMETIC_REFERENCE, $ARITHMETIC_PIXEL) evals=100
SUITE["arithmetic"]["jacobian/callable"] = @benchmarkable transform_jacobian($ARITHMETIC_CALLABLE, $ARITHMETIC_PIXEL) evals=100

SUITE["fixed_inputs"]["forward/structured"] = @benchmarkable $FIXED_REFERENCE($FIXED_PIXEL) evals=100
SUITE["fixed_inputs"]["forward/closure"] = @benchmarkable $NATIVE_CALLABLE_REFERENCE(SVector($FIXED_PIXEL[1], 300.0)) evals=100
SUITE["fixed_inputs"]["inverse/checked"] = @benchmarkable $FIXED_INVERSE($FIXED_WORLD) evals=100
SUITE["fixed_inputs"]["constant_axis"] = @benchmarkable $CONSTANT_AXIS($PIXEL) evals=100
SUITE["axis_order"]["forward/identity"] = @benchmarkable GWCS.pixel_to_world($ORDERED_3D_WCS, $PIXEL_3D) evals=100
SUITE["axis_order"]["forward/permuted"] = @benchmarkable GWCS.pixel_to_world($PERMUTED_3D_WCS, $PIXEL_3D) evals=100
SUITE["axis_order"]["inverse/identity"] = @benchmarkable GWCS.world_to_pixel($ORDERED_3D_WCS, $ORDERED_WORLD_3D) evals=100
SUITE["axis_order"]["inverse/permuted"] = @benchmarkable GWCS.world_to_pixel($PERMUTED_3D_WCS, $PERMUTED_WORLD_3D) evals=100
SUITE["domains"]["simple/inside"] = @benchmarkable GWCS.pixel_to_world($SIMPLE_DOMAIN_WCS, $DOMAIN_INSIDE) evals=100
SUITE["domains"]["simple/outside"] = @benchmarkable GWCS.pixel_to_world($SIMPLE_DOMAIN_WCS, $DOMAIN_OUTSIDE) evals=100
SUITE["domains"]["compound/inside"] = @benchmarkable GWCS.pixel_to_world($COMPOUND_DOMAIN_WCS, $DOMAIN_INSIDE) evals=100
SUITE["domains"]["compound/outside"] = @benchmarkable GWCS.pixel_to_world($COMPOUND_DOMAIN_WCS, $DOMAIN_OUTSIDE) evals=100
SUITE["domains"]["mask/simple/inside"] = @benchmarkable GWCS.in_domain($SIMPLE_DOMAIN_WCS, $DOMAIN_INSIDE) evals=100
SUITE["domains"]["mask/simple/outside"] = @benchmarkable GWCS.in_domain($SIMPLE_DOMAIN_WCS, $DOMAIN_OUTSIDE) evals=100
SUITE["domains"]["mask/compound/inside"] = @benchmarkable GWCS.in_domain($COMPOUND_DOMAIN_WCS, $DOMAIN_INSIDE) evals=100
SUITE["domains"]["mask/compound/outside"] = @benchmarkable GWCS.in_domain($COMPOUND_DOMAIN_WCS, $DOMAIN_OUTSIDE) evals=100
SUITE["domains"]["mask/image/inside"] = @benchmarkable GWCS.in_image($COMPOUND_DOMAIN_WCS, $DOMAIN_INSIDE) evals=100
SUITE["domains"]["mask/image/outside"] = @benchmarkable GWCS.in_image($COMPOUND_DOMAIN_WCS, $DOMAIN_OUTSIDE) evals=100
SUITE["domains"]["mask/batch"] = @benchmarkable GWCS.in_domain($COMPOUND_DOMAIN_WCS, $DOMAIN_BATCH) evals=1
SUITE["domains"]["grid/simple"] = @benchmarkable grid_from_bounding_box($SIMPLE_DOMAIN_BOX; fixed=(3 => 1.0,)) evals=1
SUITE["domains"]["grid/compound"] = @benchmarkable grid_from_bounding_box($COMPOUND_DOMAIN_BOX; selector=1) evals=1

if abspath(PROGRAM_FILE) == @__FILE__
    for name in keys(SUITE)
        println("\nBenchmark suite: $name")
        show_benchmarks(run(SUITE[name]; verbose=true))
    end
end
