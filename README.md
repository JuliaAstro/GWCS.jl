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
