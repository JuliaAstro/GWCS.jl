# GWCS.jl

GWCS.jl implements composable generalized world coordinate systems in Julia. Pipelines are normalized and compiled during construction so scalar evaluation uses concrete, statically sized transforms without runtime graph traversal.

## Constructing a WCS

Transforms compose left to right with `|`. Independent transforms combine with `&`.

```julia
using GWCS

detector = Frame2D(name="detector", axes_names=("x", "y"), unit=("pixel", "pixel"))
world = Frame2D(name="world", axes_names=("u", "v"))

shift = Shift(-10.0) & Shift(-20.0)
linear = AffineTransform([0.1 0.01; -0.02 0.2], [1.0, 2.0])
wcs = WCS(Step(detector, shift | linear), Step(world))

pixel_to_world(wcs, 11.0, 22.0)
world_to_pixel(wcs, 1.12, 2.38)
```

Coordinates use each boundary frame's native order. `axes_order[i]` gives the one-based transform position of native frame axis `i`; GWCS compiles static boundary remaps when this order is not the identity. Transform kernels and stored pipeline steps remain in transform order. GWCS.jl does not apply an implicit zero- or one-based pixel offset, so a FITS-equivalent pipeline must encode the FITS reference-pixel convention explicitly.

`pixel_to_world` and `world_to_pixel` accept scalar varargs, tuples, vectors, and matrices whose columns are coordinates. Bounding boxes are enforced by default when present; pass `with_bounding_box=false` to evaluate outside the box or set `fill_value` to change the default `NaN`.

### Editing a pipeline

`WCS` values are immutable. Editing functions return a newly validated and compiled WCS, leaving the source unchanged:

```julia
recalibrated = set_transform(wcs, "detector", "world", new_transform)
with_extra_step = insert_transform(recalibrated, "world", correction)
with_frame = insert_frame(with_extra_step, "detector", IdentityTransform{2}(), Frame2D(name="corrected"))
renamed = rename_frame(with_frame, "corrected", "focal")
```

`set_transform` requires adjacent frames. `insert_transform` composes before the named frame by default or after it with `after=true`; it changes coordinate values without changing frame metadata. `insert_frame(input, transform, output)` requires exactly one endpoint to already exist and injects a new edge, so it does not preserve the overall mapping unless the added transform is an identity. `remove_frame` removes an intermediate frame while composing its adjacent transforms. `prepend_transform`, `append_transform`, `prepend_frame`, and `append_frame` provide boundary conveniences. Use `get_frame(wcs, name)` to retrieve the canonical frame stored in the pipeline.

For several edits, use `WCSBuilder`. Every bang operation constructs a complete candidate WCS before committing it, so a failed edit leaves the builder unchanged:

```julia
builder = WCSBuilder(wcs)
set_transform!(builder, "detector", "world", new_transform)
insert_frame!(builder, "detector", IdentityTransform{2}(), Frame2D(name="focal"))
edited = WCS(builder)
```

The builder is a construction-time object whose stored WCS type may change after each edit. Extract the immutable result with `WCS(builder)` before performance-critical evaluation. Edits that keep the input frame preserve `pixel_shape` and `bounding_box`. Prepending a new input frame requires explicit replacements, or explicit `nothing` values to clear existing input-domain metadata because GWCS cannot infer that metadata in the new coordinate system.

Use `in_domain(wcs, pixels...)` when only the validity decision is needed. It returns `true` for finite input coordinates accepted by the WCS domain; the matrix form returns one Boolean per coordinate column. `in_image(wcs, world...)` performs the corresponding world-side test by requiring a successful finite inverse in the input domain:

```julia
in_domain(wcs, 11.0, 22.0)
in_image(wcs, pixel_to_world(wcs, 11.0, 22.0))
mask = in_domain(wcs, pixel_matrix)
```

`in_image(...; with_bounding_box=false)` tests finite invertibility without enforcing the input domain. Numerical-inverse controls such as `tolerance`, `maxiter`, `initial_guess`, and `max_backtracks` are supported. Nonconvergence returns `false`; invalid solver arguments and ambiguous multi-entry compound domains still raise errors.

### Bounding boxes and input domains

`BoundingBox` constrains all input axes by default. Use the dimensioned constructor to leave selected axes unconstrained:

```julia
box = BoundingBox{3}(((0.0, 4095.0), (0.0, 4095.0)), (1, 2))
wcs = WCS(Step(detector, transform), Step(world); bounding_box=box)

bounded_axes(box) # (1, 2)
ignored_axes(box) # (3,)
```

`CompoundBoundingBox` selects a simple box from exact values of one or more input axes. This represents detector regions, slices, and orders whose valid spatial bounds differ:

```julia
boxes = Dict(
    (1,) => BoundingBox{3}(((0.0, 2047.0), (0.0, 2047.0)), (1, 2)),
    (2,) => BoundingBox{3}(((2048.0, 4095.0), (0.0, 2047.0)), (1, 2))
)
domain = CompoundBoundingBox(boxes, (3,))
```

Selector matching is exact. An unknown selector is outside the domain and therefore produces `fill_value` during WCS evaluation. `select_bounding_box(domain, selector)` or `domain[selector]` returns one selected box explicitly and throws `KeyError` when it is absent.

Use `with_bounding_box(wcs, box)` to return a newly validated WCS with a replacement domain, or pass `nothing` to remove it. `with_pixel_shape(wcs, shape)` provides the corresponding persistent pixel-shape update. `fix_inputs(wcs, ...)` reduces simple and compound domains, including selecting and collapsing a compound box when selector inputs are fixed.

`grid_from_bounding_box` returns a coordinate matrix whose columns can be passed directly to `pixel_to_world`. Scalar or per-axis steps are supported; `center=true` uses pixel-center rounding. Compound boxes require `selector`, and other ignored axes require fixed values:

```julia
pixels = grid_from_bounding_box(wcs; selector=1, fixed=("order" => 2), step=4)
world = pixel_to_world(wcs, pixels)
```

Numerical inversion can derive an initial guess from a simple box or a compound box with one selector entry. A multi-entry compound box requires an explicit `initial_guess` because the world coordinate does not generically identify a discontinuous selector branch before solving.

These predicates test the input domain, not an axis-aligned approximation of the output footprint. Output-footprint masking remains coupled to future footprint sampling and celestial longitude/pole handling.

### Coordinate frames

Frames classify axes with the standard `AxisType` values `SPATIAL`, `SPECTRAL`, `TIME`, `STOKES`, `PIXEL`, and `UNDEFINED`; generic frames may also use custom strings. `SpectralFrame`, `TemporalFrame`, and `StokesFrame` provide one-dimensional defaults, while `CompositeFrame` combines heterogeneous outputs without changing the compiled numeric kernel:

```julia
sky = CelestialFrame(name="sky", axes_order=(2, 4))
wave = SpectralFrame(name="wave", axes_order=1, axes_names="lambda", unit="nm")
time = TemporalFrame(name="time", axes_order=3, unit="s", reference_epoch=2460000.5, time_format="jd")
world = CompositeFrame((sky, wave, time); name="world")
```

Constituent `axes_order` values are global transform positions and must form a complete permutation across the composite. Frame metadata and public WCS values remain in frame-native order. Use `axis_names`, `axis_physical_types`, `axis_units`, and `axis_types` for frame metadata, or the corresponding `pixel_axis_*` and `world_axis_*` WCS accessors.

`SpectralFrame` infers common wavelength, frequency, energy, and velocity physical types from unit strings. `TemporalFrame` preserves reference epoch, scale, format, and location metadata. `offset_to_time` and `time_to_offset` handle numeric JD/MJD origins for common time units and delegate ordinary arithmetic to compatible external time values; other numeric time formats require an optional time integration. Rich units, celestial coordinates, and astronomy time objects remain optional extension concerns.

## Transforms and inverses

The initial transform set includes identity, shift, scale, affine, constant-output, fixed-input, and axis-remapping transforms, one- and two-dimensional polynomials, spherical/cartesian conversions, precomputed 3D rotations, and all 28 FITS/WCSLIB spherical projections. `compose` and `parallel` are the function forms of `|` and `&`.

Projection algorithms are separate from transform direction. `IntermediateToNative(projection)` converts projection-plane coordinates to native spherical longitude and latitude, while `NativeToIntermediate(projection)` provides the inverse direction. Both sides use degrees:

```julia
projection = ZenithalPerspective(1.5, 15.0)
pix2sky = IntermediateToNative(projection)
sky2pix = inverse_transform(pix2sky)
```

Parameterized projections use the ASDF parameter names corresponding to the FITS projection codes. Invalid parameters raise `ArgumentError`; coordinates on singular, overlapping, or non-visible projection branches raise `DomainError`. The projection types are also shared by native FITS parsing, so an algorithm has one numeric implementation regardless of its metadata source.

`celestial_transform` combines a projection with the Paper II native-to-celestial rotation derived from a reference longitude and latitude. It computes the default `LONPOLE` when omitted and uses `LATPOLE` to disambiguate general pole solutions:

```julia
using StaticArrays

sky = celestial_transform(GnomonicProjection(), (83.8221, -5.3911))
sky(SVector(0.0, 0.0)) # reference longitude and latitude
```

This constructor returns a transform rather than a complete WCS or frame. It is the low-level construction surface used by the FITS parser and programmatic pipelines.

Compatible transforms can also be combined elementwise with `+`, `-`, `*`, `/`, and `^`. Both operands receive the same coordinates and must have matching input and output dimensions:

```julia
correction = Shift(1.0) * Scale(2.0)
correction(SVector(3.0)) # [24.0]
```

`elementwise(operator, left, right)` is the function form. Arithmetic transforms recursively compile their operands, propagate axis dependencies and Jacobians, and remain allocation-free for scalar static-vector evaluation. They do not infer an analytical inverse because combining invertible operands does not generally produce a uniquely invertible function; square WCS pipelines can use numerical inversion or an explicit `WithInverse`.

Analytical and metadata-provided inverses are compiled into the WCS. Square pipelines without an inverse use a damped Newton solver with ForwardDiff Jacobians:

```julia
world_to_pixel(wcs, world_coordinate; tolerance=1e-8, maxiter=50, initial_guess=(10.0, 20.0), quiet=false)
```

The default `quiet=true` returns the best solution on nonconvergence. `quiet=false` throws `NoConvergence`, which records the best solution, accuracy, and iteration count. Numerical inversion is not available for pipelines with different pixel and world dimensions.

### Callable transforms

`CallableTransform` is a low-boilerplate entry point for programmatic transforms. Declare the input and output dimensions, then provide functions that accept positional scalar coordinates:

```julia
using StaticArrays

forward(x, y) = (x + y, x - y)
backward(u, v) = ((u + v) / 2, (u - v) / 2)
jacobian(x, y) = SMatrix{2,2}(1, 1, 1, -1)

dependencies = SMatrix{2,2,Bool}(true, true, true, true)
transform = CallableTransform{2,2}(forward; inverse=backward, jacobian, dependencies)
wcs = WCS(Step(detector, transform), Step(world))
```

A one-output function may return a scalar. Multi-output functions return a tuple or vector of the declared length, and Jacobians return an output-by-input matrix. Optional `dependencies` and `inverse_dependencies` matrices declare structural axis correlations; omitted declarations conservatively assume every output depends on every input. The callable, inverse, Jacobian, and dependency declarations are stored in concrete type parameters, so named functions and type-stable closures can remain inferred and allocation-free. Supplied Jacobians are propagated through composed and parallel transforms and used by numerical inversion; otherwise GWCS.jl uses ForwardDiff.

Use `CallableTransform` for local functions, prototypes, and closures over calibration values. Define a dedicated `AbstractTransform` subtype when a transform needs stable metadata tags, serialization, specialized compilation, parameter introspection, or a reusable package API. Arbitrary functions are not automatically serializable, and allocations or type instability inside user callables remain visible during evaluation.

### Region-selected transforms

`RegionSelector` implements discontinuous transforms such as IFU slices, detector regions, and spectral orders. A label mapper chooses a label from the full input coordinate, then the selector evaluates the transform registered for that label:

```julia
labels = [0 2; 1 0]
mapper = ArrayLabelMapper(labels; no_label=0)
regions = RegionSelector(mapper, (
    1 => (Shift(1.0) & Scale(1.0)),
    2 => (Shift(2.0) & Scale(2.0))
))

regions(SVector(0.0, 1.0)) # [1.0, 1.0]
regions(SVector(1.0, 0.0)) # [3.0, 0.0]
```

Array lookup follows the GWCS pixel convention: `(x, y)` coordinates are zero-origin pixel centers, converted with `floor(value + 0.5)`, and used to index the label array as `[y, x]`. Nonfinite and out-of-array coordinates have no label. `DictionaryLabelMapper` supports tolerance-based numeric keys, `RangeLabelMapper` uses validated nonoverlapping open intervals, `TransformLabelMapper` adapts a one-output numeric transform, and `CallableLabelMapper{N,L}` accepts arbitrary label logic returning `L` or `nothing`.

GWCS.jl represents no match as `nothing`, so numeric zero and the empty string remain available as ordinary labels. A missing label or a label without a registered branch returns `undefined_value` in every output, defaulting to `NaN`. `selected_transform(regions, label)` retrieves one branch explicitly.

`WithInverseLabelMapper` pairs forward- and world-side label mappers. A `RegionSelector` has an analytical inverse only when the paired inverse mapper and every selected transform are invertible. Dependencies conservatively include both branch calculations and inputs used to choose the branch; Jacobians use the selected branch and raise `DomainError` where no branch applies.

Region selection and compound bounding boxes remain separate. The compound box determines whether an input is in the declared WCS domain, while the label mapper determines which transform applies there. This allows `in_domain` to report domain membership independently of whether a region transform is defined.

### Fixed inputs, constants, and axis names

A closure is sufficient for one-off partial application, but `fix_inputs` creates a lower-dimensional transform that remains visible to composition, Jacobian propagation, metadata parsing, and inversion:

```julia
transform = AnnotatedTransform(AffineTransform([1.0 1.0; 1.0 -1.0]); inputs=("x", "y"), outputs=("u", "v"))
slice = fix_inputs(transform, "y" => 2.0)

slice(SVector(3.0))                    # [5.0, 1.0]
inverse_transform(slice)(SVector(5.0, 1.0))
```

Programmatic integer indices are one-based. String and symbol keys use `input_names(transform)`. `annotate_transform` attaches an optional transform name and input/output names without changing the numerical transform:

```julia
transform = annotate_transform(transform; name="sum_and_difference")
transform_name(transform)
input_names(transform)
output_names(transform)
```

When the original transform has an inverse, the fixed-input inverse verifies both the recovered fixed coordinates and the complete forward round trip. It throws `DomainError` for an output outside the fixed slice instead of returning a misleading coordinate. `atol` and `rtol` keywords on `fix_inputs` control this validation.

`fix_inputs(wcs, ...)` also projects the input frame, bounding box, and pixel shape. Fixing every input produces a zero-dimensional WCS evaluated with `pixel_to_world(wcs)`:

```julia
detector_slice = fix_inputs(wcs, "y" => 20.0)
constant_wcs = fix_inputs(wcs, ("x" => 11.0, "y" => 20.0))
pixel_to_world(constant_wcs)
```

`ConstantTransform(value; n_inputs=0)` inserts a constant output axis without consuming an input. Constants with ignored inputs have no inverse. A zero-input constant has a checked inverse that accepts only its constant output and returns an empty coordinate.

### Axis dependencies

`axis_correlation_matrix(wcs)` returns a static boolean matrix with world axes on rows and pixel axes on columns, both in their frames' native order. `axis_correlation_matrix(transform)` reports the same output-by-input relationship for a transform:

```julia
axis_correlation_matrix(wcs)
transform_dependencies(transform)
```

Composition uses boolean dependency propagation, parallel transforms produce block-diagonal correlations, and fixed inputs remove the corresponding columns. Unknown external transforms use a conservative dense matrix. External transform packages may specialize `transform_dependencies`; a false entry must guarantee independence for every valid coordinate and parameter value rather than describe a Jacobian that happens to be zero at one point.

## FITS WCS

`parse_fits_wcs(header; key=' ', preserve_units=false, minerr=0.0)` parses one primary or alternate FITS image-header representation into an ordinary compiled `WCS` pipeline. `find_all_fits_wcs(header; ...)` returns every representation as a `Dict{Char,WCS}` and reports malformed alternates with their key. Parsing is separate from tagged GWCS/ASDF metadata, which continues to use `parse_wcs`.

```julia
header = Dict(
    "NAXIS" => 2, "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN",
    "CRPIX1" => 512.0, "CRPIX2" => 512.0,
    "CRVAL1" => 83.8221, "CRVAL2" => -5.3911,
    "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4,
)
wcs = parse_fits_wcs(header)
pixel_to_world(wcs, 512.0, 512.0)
```

FITS 1-based pixel coordinates are represented explicitly by the `CRPIX` affine stage; no origin offset is applied elsewhere. Public outputs retain FITS axis order, including interleaved celestial, spectral, time, Stokes, and linear axes. Canonical outputs use degrees for celestial axes, SI units for supported linear spectral axes, and seconds for time axes. Set `preserve_units=true` to retain supported `CUNIT` values at the output boundary. Incomplete image shapes use `nothing` for WCS axes without a corresponding `NAXISn`.

The parser supports PC, CD, and legacy CROTA linear forms, all 28 FITS spherical projections and their PV parameters, pole parameters, SIP, TPV, pre-2012 SCAMP TAN-plus-PV headers, Paper IV prior and sequent TPD polynomials, legacy DSS and WAT-encoded TNX/ZPX polynomials normalized to TPD at construction, Paper IV D2IM/CPDIS lookup arrays, multilinear `-TAB` coordinate arrays, Paper III spectral algorithms, and GRI/GRA grism models. SIP AP/BP and Paper IV TPD.REV coefficients are approximate reverse initializers only: `world_to_pixel` seeds the whole-WCS numerical inverse with them and always refines against the forward WCS. FITSIO and FITSFiles are weak dependencies that add detached-header and HDU adapters; linked HDU inputs copy referenced images and table columns into backend-neutral transforms during construction. Separable monotonic 1D tables have direct inverses, while coupled tables use the checked WCS numerical inverse. RESTFRQ/RESTWAV and spectral reference-system metadata are retained on spectral frames, and time reference metadata is retained on temporal frames.

## Metadata

`parse_wcs` accepts a tagged YAML string or a tagged metadata mapping. The core parser supports:

- GWCS WCS, step, simple and compound bounding boxes, generic, 2D, celestial, spectral, temporal, Stokes, and composite frames
- ICRS, structured temporal reference metadata, units, axis classifications, and nonidentity frame-axis ordering
- compose, concatenate, elementwise arithmetic, identity, shift, scale, affine, axis remapping, and polynomial transforms
- legacy and current constant transforms, fixed-input compounds, and transform names/input/output names
- array-, dictionary-, range-, and transform-backed label mappers and region selectors
- TAN, AZP, SZP, STG, SIN, ARC, ZEA, and AIR projections, native/celestial rotations, Cartesian and spherical rotation sequences, and spherical/cartesian conversions
- serialized explicit inverses

ASDF.jl is an optional dependency. Loading ASDF activates the package extension that adapts tagged nodes and lazy coefficient arrays:

```julia
using ASDF
using GWCS

file = ASDF.load("roman_cal.asdf"; extensions=true, validate_checksum=false)
wcs = parse_wcs(file["roman"]["meta"]["wcs"])
```

Unknown tags, unsupported schema versions, malformed dimensions, and missing required properties raise path-aware errors. Built-in parsers register the exact tag versions present in the corresponding ASDF schema packages instead of assuming that every future `1.x` version is compatible. Version-specific validation runs before construction, and legacy representations such as pre-1.4 constants are normalized to the current internal form. Unit-bearing transform parameters require a future units extension and are not silently stripped in the numeric core.

### Serializing metadata

`to_metadata` returns a canonical tagged tree for supported transforms, frames, bounding boxes, steps, WCS pipelines, label mappers, and region selectors. The result can be inspected directly, rendered with YAML.jl, or parsed back through the same validated readers:

```julia
using YAML

tree = to_metadata(wcs)
text = YAML.write(tree)
reparsed = parse_wcs(tree)
```

Serialization represents the compiled Julia object, not the exact source tree. Construction-time fusion may therefore replace a sequence with an equivalent affine transform, and output uses the latest exact schema versions supported by the writer. Projection trees preserve their direction and canonical numeric parameters. Unknown optional source properties and original tag versions are not preserved yet. Unsupported values such as callable closures, callable label mappers, nonstandard affine dimensions, and unregistered external transforms raise explicit errors instead of producing incomplete metadata.

Loading ASDF.jl adds `write_asdf`. It converts canonical tagged nodes to ASDF.jl's container representation and delegates file and block handling to ASDF.jl:

```julia
using ASDF

write_asdf("calibration.asdf", wcs)
write_asdf("product.asdf", Dict("meta" => Dict("wcs" => wcs), "data" => image))
```

The WCS form writes under the key `"wcs"` by default; use `key` to choose another key. Array-valued transform and selector metadata is currently written inline. Lazy, external, compressed, and explicit block-placement policies are deferred to a later ASDF integration.

### Registering external transforms

External packages can add metadata support for native `AbstractTransform` subtypes without modifying GWCS.jl. Dimensions are encoded in the subtype, while analytical inverses, Jacobians, and optimized execution forms use ordinary Julia method extensions:

```julia
module ExampleTransforms

using GWCS
using StaticArrays
import GWCS: compile_transform, inverse_transform, transform_dependencies, transform_jacobian

struct OffsetTransform{T<:Real} <: AbstractTransform{1,1}
    offset::T
end

@inline (transform::OffsetTransform)(coords::StaticVector{1}) = SVector(coords[1] + transform.offset)
inverse_transform(transform::OffsetTransform) = OffsetTransform(-transform.offset)
transform_jacobian(transform::OffsetTransform, coords::StaticVector{1}) = SMatrix{1,1}(one(eltype(coords)))
transform_dependencies(::OffsetTransform) = SMatrix{1,1,Bool}(true)
compile_transform(transform::OffsetTransform) = Shift(transform.offset)

function parse_offset(node, context)
    offset = metadata_value(required_metadata(node, "offset", context))
    return OffsetTransform(offset)
end

serialize_offset(transform, context) = Dict("offset" => transform.offset)

function __init__()
    register_transform!("tag:example.org:gwcs/offset", parse_offset; versions=(v"1.0.0", v"1.1.0"))
    register_transform_serializer!(OffsetTransform, "tag:example.org:gwcs/offset", serialize_offset; version=v"1.1.0")
end

end
```

Tags use their full namespace and omit the trailing schema version. Prefer `versions` when support is known for exact schema tags. Half-open `min_version`/`max_version` intervals remain available when an external schema explicitly guarantees compatibility, and overlapping registrations are rejected. Optional `validate` and `normalize` callbacks receive the source node and `TransformParseContext`; validation runs first, then normalization may return a canonical mapping for the constructor without changing the original metadata. A serializer receives the transform and a write context and returns an untagged mapping; its exact `version` must have a matching parser registration. Register parser and serializer support from an extension package's `__init__` so the process-wide default registry is populated at runtime. For deterministic or restricted conversion, use `registry = copy(default_transform_registry())`, register additional transforms into that isolated registry, and pass it to `parse_transform`, `parse_wcs`, or `to_metadata`.

Registration and metadata callbacks run only while constructing a transform. A parsed `WCS` stores concrete compiled transforms, so scalar evaluation performs no registry lookup or dynamic metadata dispatch. Registered constructors must return an `AbstractTransform`; custom inverse and compilation methods are checked for compatible dimensions before the pipeline is accepted.

## Tests, benchmarks, and validation

Run the package tests with:

```sh
julia --project -e 'using Pkg; Pkg.test()'
```

The benchmark suite compares equivalent affine, TAN, TAN-SIP, nonlinear spectral, grism, and spherical-projection scalar transforms against FITSWCS.jl. It also covers direct lookup-table and `-TAB` evaluation:

```sh
julia --project=benchmark benchmark/benchmarks.jl
```

External Roman validation requires the local ASDF.jl checkout configured in `validation/Project.toml`, Python GWCS in `~/.venvs/311_main`, and one or more Roman L2 ASDF files:

```sh
julia --project=validation validation/validate_asdf_parsers.jl
julia --project=validation validation/compare_asdf_parsers.jl
julia --project=validation validation/compare_arithmetic.jl
julia --project=validation validation/compare_bounding_boxes.jl
julia --project=validation validation/validate_asdf_bounding_boxes.jl
julia --project=validation validation/check_schema_versions.jl
julia --project=validation validation/compare_roman.jl FILE...
```

The parser validation checks ASDF.jl lazy-array orientation and compares identity, affine, TAN, native/celestial rotation, Euler rotation, spherical rotation-sequence, arithmetic compound, bounding-box, label-mapper, and region-selector behavior against Astropy and Python GWCS. `check_schema_versions.jl` verifies that the built-in exact-version table matches the installed ASDF schema packages. Roman validation compares four detector positions per file in both directions against Python GWCS. Its current acceptance limits are `1e-10` degrees forward and `1e-8` pixels inverse.
