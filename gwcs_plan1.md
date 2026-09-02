# GWCS.jl implementation plan

## Problem and scope

Build the first usable architecture for a pure-Julia generalized world coordinate
system package. The initial release should make transformation composition easy,
compile user- or metadata-defined pipelines into fast scalar kernels, expose
`pixel_to_world` and `world_to_pixel`, parse the transform subset used by Roman L2
ASDF metadata, and validate correctness and performance against FITSWCS.jl and
Python GWCS.

The target repository currently contains only `README.md`, `.gitignore`, and
license/git metadata, so package structure and APIs will be established from
scratch.

## Reference-code findings

- FITSWCS.jl uses `StaticArrays`, concretely typed WCS fields, tuple/vararg scalar
  entry points, and `SVector` intermediates. Its benchmark project uses
  `BenchmarkTools` suites grouped by forward, inverse, and construction work.
- Python GWCS represents a WCS as ordered frame/transform `Step` objects. It
  dynamically combines Astropy models for forward or backward traversal.
- Astropy model composition has two essential operators:
  sequential composition (`|`) and independent/parallel concatenation (`&`).
  Explicit user inverses are serialized in an `inverse` property.
- Representative Roman L2 files all contain four frames:
  `detector -> v2v3 -> v2v3vacorr -> world`.
- The observed Roman WCS subset uses these tags:
  `wcs`, `step`, `frame2d`, `celestial_frame`, ICRS, unit, bounding box,
  `compose`, `concatenate`, `shift`, `scale`, `remap_axes`, `polynomial`,
  `rotate_sequence_3d`, and `spherical_cartesian`.
- Roman detector models contain explicit inverse subtrees. The implementation
  must preserve these instead of falling back to iterative inversion.
- The local ASDF.jl checkout can read the files with
  `ASDF.load(path; extensions=true, validate_checksum=false)` and retains unknown
  tags in `TaggedMapping`, `TaggedSequence`, and `TaggedScalar` nodes.

## Architectural decisions

### 1. Separate construction from execution

Use a two-layer design:

1. A construction layer accepts nested transform objects, pipeline steps, and
   tagged metadata. It may use dynamic registries and spend additional time
   validating, normalizing, and compiling the graph.
2. An immutable execution layer stores only concretely typed transforms and
   precompiled forward/inverse chains. `pixel_to_world` and `world_to_pixel`
   perform no frame lookup, graph traversal, registry access, or composition.

This permits ergonomic construction without putting dynamic dispatch in the hot
path.

### 2. Core transform interface

Define:

```julia
abstract type AbstractTransform{Nin,Nout} end

input_dimension(::AbstractTransform{Nin}) where {Nin} = Nin
output_dimension(::AbstractTransform{Nin,Nout}) where {Nin,Nout} = Nout
```

Every compiled transform is callable on `SVector{Nin,T}` and returns an
`SVector{Nout,T′}`. Leaf methods and composition evaluators are `@inline`.
Abstract-vector, tuple, and scalar-vararg APIs convert once at the outer function
barrier.

Initial leaf types:

- `IdentityTransform`
- `Shift`
- `Scale`
- `AffineTransform` with static matrix/vector storage
- `RemapAxes`
- `Polynomial1D` and `Polynomial2D`, with degree encoded in the type and
  coefficients stored statically
- `SphericalToCartesian` and `CartesianToSpherical`
- `Rotation3D`, storing one precomputed `SMatrix{3,3}` rather than evaluating a
  sequence of rotations on every call

Initial combinators:

- `ComposedTransform{Ts<:Tuple}` for flattened sequential transforms
- `ParallelTransform{Ts<:Tuple}` for independent transforms over disjoint input
  slices
- `WithInverse{F,I}` for metadata-provided inverse pairs

Expose `compose(...)` and `parallel(...)` as the canonical constructors. Add
transform-specific `|` and `&` convenience methods to mirror Astropy metadata and
make programmatic translation straightforward.

### 3. Compilation and optimization

`compile_transform` will:

- validate dimensions at every edge;
- flatten nested sequential and parallel nodes;
- remove identities;
- fuse adjacent scalar shifts/scales and parallel affine operations;
- precompute rotation-sequence matrices and their transposes;
- convert polynomial coefficient arrays to static storage and preserve Astropy's
  `c[i,j] * x^i * y^j` orientation;
- compile forward and inverse roots once;
- retain explicit ASDF inverse subtrees through `WithInverse`;
- attach a typed `NoInverse` marker when an inverse cannot be compiled, allowing
  `world_to_pixel` to select the numerical fallback without exception-driven
  control flow.

Analytical and metadata-provided inverses remain the fast path. For a square
pipeline without one, `world_to_pixel` will automatically invoke a numerical
solver. The first implementation will use `ForwardDiff` to form the Jacobian of
the compiled scalar kernel, then apply damped Newton iterations with
backtracking, finite-value/divergence checks, wrapped residuals for longitude
axes, and a bounding-box-center default initial guess. Expose tolerance,
iteration limit, initial guess, and quiet/error behavior as keywords. Match
Python GWCS by defaulting to `quiet=true`, which returns the best solution;
`quiet=false` raises a `NoConvergence` exception carrying the best solution,
residual/step accuracy, and iteration count. Non-square transforms will fail
explicitly.

### 4. Frames, steps, and WCS pipeline

Frame metadata stays outside the numerical kernel:

- `AbstractCoordinateFrame{N}`
- `CoordinateFrame` / `Frame2D`
- `CelestialFrame`
- fields for name, axes names/order, physical types, unit strings, and reference
  frame metadata

Pipeline construction uses:

```julia
Step(frame, transform)
Step(output_frame) # final step, no transform
WCS(steps...; name="", pixel_shape=nothing, bounding_box=nothing)
```

Construction validates unique frame names, matching dimensions, a terminal
transform-free step, and invertibility metadata. The resulting `WCS` stores:

- immutable steps for inspection;
- a compiled full forward transform;
- a compiled full reverse transform or typed `NoInverse` marker;
- optional typed bounding-box metadata;
- input/output frame metadata.

Provide a slower, explicit `transform(wcs, from_frame, to_frame, coords...)` for
intermediate-frame use. It may construct and return a compiled subtransform;
the main forward/reverse calls never use this path.

Pixel coordinates are pipeline-native. No hidden zero/one-based correction is
performed by `pixel_to_world` or `world_to_pixel`. Roman metadata therefore
matches Python GWCS's detector coordinates, while a programmatic FITS-equivalent
pipeline explicitly encodes FITS's 1-based reference-pixel shift.

### 5. Public evaluation API

Match FITSWCS.jl call shapes:

```julia
pixel_to_world(wcs, coords::AbstractVector)
pixel_to_world(wcs, coords::Tuple)
pixel_to_world(wcs, x::Real, y::Real, ...)
world_to_pixel(wcs, coords::AbstractVector)
world_to_pixel(wcs, coords::Tuple)
world_to_pixel(wcs, lon::Real, lat::Real, ...)
```

Add matrix batch methods where each column is one coordinate, implemented as
tight loops over the scalar compiled kernel. Preserve `Float32` where all model
parameters can be converted safely; use promoted floating-point output
otherwise.

Support `with_bounding_box` and `fill_value` keywords. A WCS without a bounding
box should compile to a no-op path. Defaults should match imported Python GWCS
behavior, while programmatically constructed WCS objects without a box pay no
runtime cost. `world_to_pixel` additionally accepts numerical-inverse controls
such as `tolerance`, `maxiter`, `initial_guess`, and `quiet`; they are ignored
when a compiled inverse exists.

### 6. Metadata parser and ASDF integration

Use `StaticArrays`, `YAML`, and `ForwardDiff` as direct dependencies, plus stdlib
`LinearAlgebra`. Keep ASDF.jl as a weak dependency with a package extension.

Core parsing API:

```julia
parse_wcs(metadata::AbstractString)
parse_wcs(metadata::AbstractDict)
parse_transform(metadata)
```

The YAML string parser will preserve unknown/custom tags in small GWCS-owned
tagged-node wrappers. Parsing uses a construction-time `TagRegistry` keyed by
normalized tag family and schema version. Unknown tags, malformed dimensions,
unsupported versions, and missing required properties must fail with the tag
and metadata path in the error.

Implement parsers for the exact Roman subset listed above, including:

- common transform properties: `name`, `inputs`, `outputs`, `inverse`, and
  `bounding_box`;
- composition and concatenation;
- inline or array-backed polynomial coefficients;
- frame/unit metadata;
- Cartesian rotation sequences and spherical/cartesian conversion.

`ext/GWCSASDFExt.jl` will adapt ASDF `Tagged*` values and materialize
`ASDF.NDArray` coefficients before passing a package-owned metadata tree to the
core parser. ASDF remains optional for ordinary programmatic or YAML-string use.

The normal test suite will parse an extracted Roman WCS YAML fixture with inline
coefficient arrays, so it does not depend on ASDF.jl or local science files.

### 7. Proposed repository layout

```text
Project.toml
src/
  GWCS.jl
  frames.jl
  transforms/
    core.jl
    elementary.jl
    composition.jl
    polynomial.jl
    geometry.jl
  pipeline.jl
  api.jl
  metadata.jl
ext/
  GWCSASDFExt.jl
test/
  runtests.jl
  test_transforms.jl
  test_composition.jl
  test_pipeline.jl
  test_api.jl
  test_metadata.jl
  fixtures/roman_wcs.yaml
benchmark/
  Project.toml
  benchmarks.jl
validation/
  python_reference.py
  compare_roman.jl
```

## Implementation todos

1. **Scaffold the Julia package.**
   Create `Project.toml`, module entry point, source/test directories, Julia
   1.10 compatibility, direct dependencies, ASDF weak dependency/extension, and
   basic package-load tests.
2. **Implement typed leaf transforms and inverses.**
   Add dimensions, static scalar evaluation, affine/remap/polynomial/geometry
   models, explicit inverse handling, validation, and focused numerical tests
   against Astropy formulas.
3. **Implement composition and compilation.**
   Add sequential and parallel combinators, normalization/fusion, type-stable
   tuple recursion, inverse-chain construction, and zero-allocation scalar
   checks.
4. **Implement frames, steps, and compiled WCS pipelines.**
   Add pipeline validation, frame lookup for construction/inspection,
   intermediate-frame transforms, bounding boxes, and immutable forward/reverse
   execution roots.
5. **Implement the FITSWCS-style API.**
   Add scalar varargs, tuples, vectors, and matrix batches for
   `pixel_to_world`/`world_to_pixel`, dimension errors, floating-point
   promotion, and allocation/type-inference tests.
6. **Implement automatic numerical inversion.**
   Add ForwardDiff-compatible scalar kernels, damped Newton/backtracking for
   square pipelines, longitude-aware residuals, bounding-box initialization,
   convergence diagnostics, batch behavior, and `NoConvergence` tests.
7. **Implement tagged YAML metadata parsing.**
   Add the registry, schema-version dispatch, common transform properties, the
   Roman tag subset, contextual errors, and a copied Roman YAML fixture with
   inline coefficients.
8. **Implement optional ASDF.jl integration.**
   Add the extension adapter, local-checkout development setup, lazy-array
   materialization, and a non-default integration test/script using
   `ASDF.load(...; extensions=true, validate_checksum=false)`.
9. **Build comparative benchmarks.**
   Follow FITSWCS.jl's benchmark layout and compare equivalent 2D affine, TAN,
   and TAN plus polynomial-distortion pipelines in both directions, using
   `SVector`, tuple, and ordinary-vector scalar inputs. Report time,
   allocations, and GWCS/FITSWCS ratios; benchmark construction separately.
10. **Cross-validate Roman files against Python GWCS.**
   Use `~/.venvs/311_main/bin/python` to generate forward and inverse reference
   values for multiple detectors/filters and interior/boundary points. Load the
   same files with the local ASDF.jl checkout, parse their WCS nodes, and compare
   Julia results to Python. Record reproducible commands and tolerances.
11. **Document the architecture and supported subset.**
    Expand the README with construction examples, coordinate-origin policy,
    supported ASDF tags/schema versions, inverse behavior, optional ASDF usage,
    benchmark invocation, and known first-release limits.

## Validation and acceptance criteria

- Package loads and tests pass on Julia 1.10 and current stable Julia.
- Scalar `SVector` forward and analytical/explicit inverse calls are inferred
  and allocate zero bytes after compilation.
- Equivalent affine, TAN, and distorted-TAN pipelines agree with FITSWCS.jl
  within floating-point tolerance in both directions.
- Performance comparison is ratio-based to avoid machine-specific nanosecond
  gates. Target zero allocations and median scalar time within 25 percent of
  FITSWCS.jl for affine/TAN equivalents, with no more than 50 percent overhead
  accepted initially for the generic distorted pipeline. Simple affine
  transforms should retain a stretch target below 100 ns where the machine and
  Julia version permit it.
- The extracted Roman YAML fixture parses without ASDF.jl and evaluates both
  directions.
- A square pipeline without an analytical or explicit inverse automatically
  converges on representative affine, polynomial-distortion, and celestial
  round trips. Non-convergence returns the best solution by default and raises
  structured diagnostics with `quiet=false`; non-square pipelines fail
  explicitly.
- External validation covers at least three Roman files spanning different WFI
  detectors and filters, with several detector points per file. Target agreement
  with Python GWCS at approximately `1e-10` degrees forward and `1e-8` pixels
  backward, tightening tolerances when observed numerical behavior permits.
- Unknown tags, unsupported schema versions, dimension mismatches, and absent
  inverses produce explicit errors rather than silent fallbacks.

## Notes and risks

- Deeply nested transform types can increase compilation latency. Flattening
  chains and parallel groups into tuples limits nesting while retaining static
  dispatch. Construction latency is intentionally traded for evaluation speed.
- Degree-5 Roman polynomials need a Horner-style evaluator and static
  coefficients; direct repeated exponentiation is not acceptable for the hot
  path.
- ASDF coefficient storage is row/column sensitive. Tests must use asymmetric
  coefficients so transposition bugs cannot pass.
- Rotation order and angle sign must be cross-checked against Astropy's
  `RotationSequence3D`; precomputing a matrix is valid only after matching its
  convention exactly.
- Forward-mode automatic differentiation requires every scalar transform kernel
  to accept generic real element types. Angle wrapping and bounding-box branches
  must be kept outside Jacobian evaluation or handled with explicit residual
  logic near discontinuities.
- Timing should not be a hard CI pass/fail criterion. Inference, allocations,
  correctness, and generated benchmark reports are stable regression signals;
  absolute latency is environment-dependent.
