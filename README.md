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

Coordinates are pipeline-native. GWCS.jl does not apply an implicit zero- or one-based pixel offset. A FITS-equivalent pipeline must encode the FITS reference-pixel convention explicitly.

`pixel_to_world` and `world_to_pixel` accept scalar varargs, tuples, vectors, and matrices whose columns are coordinates. Bounding boxes are enforced by default when present; pass `with_bounding_box=false` to evaluate outside the box or set `fill_value` to change the default `NaN`.

## Transforms and inverses

The initial transform set includes identity, shift, scale, affine and axis-remapping transforms, one- and two-dimensional polynomials, spherical/cartesian conversions, precomputed 3D rotations, and TAN/native/celestial projection transforms. `compose` and `parallel` are the function forms of `|` and `&`.

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

transform = CallableTransform{2,2}(forward; inverse=backward, jacobian)
wcs = WCS(Step(detector, transform), Step(world))
```

A one-output function may return a scalar. Multi-output functions return a tuple or vector of the declared length, and Jacobians return an output-by-input matrix. The callable, inverse, and Jacobian are stored in concrete type parameters, so named functions and type-stable closures can remain inferred and allocation-free. Supplied Jacobians are propagated through composed and parallel transforms and used by numerical inversion; otherwise GWCS.jl uses ForwardDiff.

Use `CallableTransform` for local functions, prototypes, and closures over calibration values. Define a dedicated `AbstractTransform` subtype when a transform needs stable metadata tags, serialization, specialized compilation, parameter introspection, or a reusable package API. Arbitrary functions are not automatically serializable, and allocations or type instability inside user callables remain visible during evaluation.

## Metadata

`parse_wcs` accepts a tagged YAML string or a tagged metadata mapping. The core parser supports the Roman L2 subset used by the validation fixtures:

- GWCS WCS, step, bounding-box, 2D frame, celestial frame, ICRS, and unit tags
- compose, concatenate, shift, scale, axis remapping, and polynomial transforms
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
import GWCS: compile_transform, inverse_transform, transform_jacobian

struct OffsetTransform{T<:Real} <: AbstractTransform{1,1}
    offset::T
end

@inline (transform::OffsetTransform)(coords::StaticVector{1}) = SVector(coords[1] + transform.offset)
inverse_transform(transform::OffsetTransform) = OffsetTransform(-transform.offset)
transform_jacobian(transform::OffsetTransform, coords::StaticVector{1}) = SMatrix{1,1}(one(eltype(coords)))
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
