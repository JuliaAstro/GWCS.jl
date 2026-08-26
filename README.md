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

## Transforms and inverses

The initial transform set includes identity, shift, scale, affine, constant-output, fixed-input, and axis-remapping transforms, one- and two-dimensional polynomials, spherical/cartesian conversions, precomputed 3D rotations, and TAN/native/celestial projection transforms. `compose` and `parallel` are the function forms of `|` and `&`.

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

## Metadata

`parse_wcs` accepts a tagged YAML string or a tagged metadata mapping. The core parser supports the Roman L2 subset used by the validation fixtures:

- GWCS WCS, step, bounding-box, 2D frame, celestial frame, ICRS, unit, and nonidentity frame-axis-order metadata
- compose, concatenate, shift, scale, axis remapping, and polynomial transforms
- legacy and current constant transforms, fixed-input compounds, and transform names/input/output names
- 3D rotation sequences and spherical/cartesian conversions
- serialized explicit inverses

ASDF.jl is an optional dependency. Loading ASDF activates the package extension that adapts tagged nodes and lazy coefficient arrays:

```julia
using ASDF
using GWCS

file = ASDF.load("roman_cal.asdf"; extensions=true, validate_checksum=false)
wcs = parse_wcs(file["roman"]["meta"]["wcs"])
```

Unknown tags, unsupported major schema versions, malformed dimensions, and missing required properties raise explicit errors. The current parser is intentionally limited to the observed Roman transform subset rather than the full Astropy modeling schema.

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

function __init__()
    register_transform!("tag:example.org:gwcs/offset", parse_offset; min_version=v"1.0.0", max_version=v"2.0.0")
end

end
```

Tags use their full namespace and omit the trailing schema version. Version intervals are half-open, and overlapping registrations are rejected. Register from an extension package's `__init__` so the process-wide default registry is populated at runtime. For deterministic or restricted parsing, use `registry = copy(default_transform_registry())`, register additional transforms into that isolated registry, and pass it to `parse_transform` or `parse_wcs`.

Registration and metadata callbacks run only while constructing a transform. A parsed `WCS` stores concrete compiled transforms, so scalar evaluation performs no registry lookup or dynamic metadata dispatch. Registered constructors must return an `AbstractTransform`; custom inverse and compilation methods are checked for compatible dimensions before the pipeline is accepted.

## Tests, benchmarks, and validation

Run the package tests with:

```sh
julia --project -e 'using Pkg; Pkg.test()'
```

The benchmark suite compares equivalent affine, TAN, and TAN-SIP scalar transforms against FITSWCS.jl:

```sh
julia --project=benchmark benchmark/benchmarks.jl
```

External Roman validation requires the local ASDF.jl checkout configured in `validation/Project.toml`, Python GWCS in `~/.venvs/311_main`, and one or more Roman L2 ASDF files:

```sh
julia --project=validation validation/compare_roman.jl FILE...
```

The validation compares four detector positions per file in both directions against Python GWCS. The current acceptance limits are `1e-10` degrees forward and `1e-8` pixels inverse.
