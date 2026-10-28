# AGENTS.md

These instructions apply to the entire GWCS.jl repository.

## Design priorities

GWCS.jl is a compiled, composable numeric WCS implementation. Optimize for:

1. Correct forward and inverse coordinate results.
2. Inferred, allocation-free scalar evaluation.
3. Static dimensions and concrete execution types.
4. Easy programmatic construction and extension.
5. Faithful parsing of Python GWCS and ASDF metadata.

Pipeline construction may perform extra validation, normalization, fusion, precomputation, and compilation when that makes repeated evaluation faster. Keep metadata lookup, schema interpretation, frame bookkeeping, units, and dynamic graph traversal out of the scalar execution path.

## Performance invariants

- Treat `pixel_to_world` and `world_to_pixel` scalar calls as hot paths.
- Preserve type inference and zero allocations for representative analytical scalar transforms after warm-up.
- Use `StaticVector` and `StaticMatrix` values in scalar kernels.
- Encode input and output dimensions in types through `AbstractTransform{Nin,Nout}`.
- Store concrete transforms in tuples or other specialization-friendly structures. Avoid abstractly typed transform collections in evaluation code.
- Compile or canonicalize transforms during `WCS` construction rather than interpreting a mutable graph during evaluation.
- Precompute matrices, permutations, inverse data, dependency information, and other reusable values.
- Benchmark dispatch or abstraction layers against direct transform evaluation when adding them.
- Do not trade clear correctness for speculative micro-optimizations. Add generated functions or complex specialization only when inference or measurements justify them.
- Keep high-level coordinate construction, units, and metadata objects outside the compiled numeric kernel.
- Construction and first-call compile latency matter as well as per-call speed. Avoid designs that multiply pipeline types per WCS, such as storing several specialized copies of the same pipeline.

## Transform architecture

- New numeric transforms should subtype `AbstractTransform{Nin,Nout}`.
- Implement scalar evaluation as `transform(coords::StaticVector{Nin}) -> StaticVector{Nout}`.
- Extend the established protocol as appropriate:
  - `inverse_transform`
  - `compile_transform`
  - `transform_jacobian`
  - `transform_dependencies`
  - `transform_name`
  - `input_names`
  - `output_names`
- Return `NoInverse()` when no analytical or explicit inverse exists. Do not invent an inverse merely because each component appears invertible.
- Return `NaN` in every output component where a coordinate has no image (an undefined analytical inverse, a point off a projection, an output outside a constrained slice); do not throw `DomainError` from transform evaluation. The evaluation API turns finite input that yields `NaN` into `fill_value`, or into a diagnosed `DomainError` with `strict=true`. Zero-dimensional results, which cannot carry `NaN`, are the only exception.
- Outer constructors should compute type parameters and delegate inward so all call forms apply the same conversions.
- Parametric inner constructors must accept unconstrained values and let field declarations or `new` perform conversion. Do not require arguments to already have the field type.
- GWCS.jl is unreleased, so breaking a public API is acceptable when it clearly simplifies the design. Propose any API-breaking change and get the maintainer's approval before making it; never break an API on your own initiative. The README documents what the package currently does, not its history, so describe the new behavior there rather than the change.

Composition with `|` is left to right. Parallel composition with `&` operates on consecutive input slices. Elementwise arithmetic sends the same inputs to both operands. Preserve these semantics when simplifying or compiling transform trees.

## Label mappers and discontinuous transforms

Labels are dispatch metadata, not numeric coordinates. Keep `AbstractLabelMapper` separate from `AbstractTransform`.

- A label mapper returns a concrete label or `nothing`.
- `nothing` means no match; do not reserve `0` or `""` as universal Julia sentinels.
- `RegionSelector` evaluates to `NaN` when no branch applies. Its `undefined_value` (Python GWCS's `undefined_transform_value`) is stored only so ASDF files round-trip exactly, and the inverse selector does not inherit it, as in Python GWCS.
- Keep region selection separate from input domains. A `CompoundBoundingBox` decides domain membership; a label mapper decides which transform applies.
- Match Python GWCS semantics deliberately when parsing its metadata. For example, dictionary label matching uses NumPy's asymmetric `isclose` rule, not Julia's `isapprox` rule.
- Construct an inverse region selector only when the inverse label mapper and every branch provide compatible inverses.

## Frames, coordinates, and domains

Transform order is the order of a transform's inputs and outputs, and it is the order of every public coordinate. Frame order is the order in which a frame stores its own axes: a `CelestialFrame` is (lon, lat), and a `CompositeFrame` concatenates its constituents. `axes_order[i]` is the one-based transform position of the frame's axis `i`, so it labels coordinates and never reorders them; a constituent's entries are positions in the whole composite. This is Python GWCS's convention.

- Public coordinates, the axis accessors (`axis_names`, `world_axis_physical_types`, and so on), bounding boxes, `pixel_shape`, and `fix_inputs` indices use transform order. Frame fields and constructor arguments use frame order.
- ASDF per-axis frame lists are in frame order, which is how Python GWCS reads them and, since gwcs#687, writes them.
- Do not add an implicit zero- or one-based pixel offset. Encode external pixel-origin conventions explicitly at interoperability boundaries.
- Matrix evaluation uses one coordinate per column.
- `in_domain` tests finite input coordinates against the input domain.
- `in_image` requires successful finite inversion and optionally enforces the input domain.
- Ignored bounding-box axes are geometrically unconstrained but must still be finite for validity predicates.
- Unknown compound selectors are outside the domain.
- Keep input-domain masks distinct from future output-footprint geometry.

## Metadata and interoperability

- Metadata parsing is a construction-time concern. Parsed WCS objects must contain concrete compiled transforms.
- Produce path-aware errors for malformed metadata.
- Validate schema versions explicitly; do not assume all future `1.x` tags are compatible.
- Use `register_transform!` for external `AbstractTransform` metadata support.
- Do not force label mappers into the transform-registration API.
- Preserve ASDF.jl as a weak dependency and adapt its tagged or lazy values through the package extension.
- GWCS.jl owns FITS WCS keyword interpretation (`src/fits.jl`); delegate FITS container I/O to FITSIO.jl or FITSFiles.jl and ASDF container I/O to ASDF.jl. Do not reimplement those container formats here.
- When behavior is intended to match Python GWCS, Astropy, NumPy, FITSWCS.jl, or WCSLIB, verify the exact convention rather than substituting a similar Julia API.
- For FITS WCS semantics, WCSLIB is the source of truth, settled from its source or a direct WCSLIB run. Always cross-check astropy as well, since it normalizes headers and implements SIP and Paper IV distortions itself; when the two disagree, state which one GWCS.jl matches and why.
- A serialization round trip must reproduce the WCS mathematically, not the source tree. Preserve frames, annotations (names, inputs, outputs), explicit inverses, and bounding boxes, and change values at most by rounding. `compose` and `parallel` may canonicalize unannotated transforms, so every form a stored step can take needs a serializer that writes standard tags Python GWCS reads.
- Keep rich units and high-level coordinate objects in optional extensions unless the numeric core requires them.

## Error handling

- Fail explicitly on invalid dimensions, ambiguous metadata, incompatible inverses, and unsupported operations.
- Do not use broad catches, silent defaults, or success-shaped fallbacks. Returning `fill_value` for coordinates without a result is the deliberate APE 14-style default, not a fallback; fail explicitly on invalid inputs, parameters, and metadata.
- Catch only errors that have defined public meaning. For example, `in_image` converts the `DomainError` of a zero-dimensional inverse to `false` (it reads numerical non-convergence from the solver's status rather than catching `NoConvergence`), and `CallableTransform` converts a `DomainError` from user code to `NaN`.
- Preserve useful diagnostic information in errors such as `NoConvergence`.

## Julia style and formatting

FITSWCS.jl at `~/Development/julia/FITSWCS.jl` is the authoritative formatting reference.

- Use American spelling.
- Keep function calls and signatures on one line when they reasonably fit.
- Avoid gratuitous vertical argument lists, line breaks, blank lines, and trailing commas.
- Prefer compact, idiomatic Julia over deeply layered helper abstractions.
- Use short keyword syntax when the variable name matches: `f(; maxiter)`, not `f(; maxiter=maxiter)`.
- Comment only non-obvious logic. Explain why a step exists rather than restating the code.
- Follow nearby naming, docstring, export, and file-layout conventions.
- Do not run formatters or broad whitespace normalization unless explicitly requested.
- Keep diffs focused on the requested behavior.

## Dependencies and compatibility

- Julia 1.10 is the minimum supported release. New syntax and dependencies must work there.
- Do not add dependencies unless they clearly reduce maintenance burden or are explicitly required.
- When adding a package, add a bounded `[compat]` entry compatible with Julia 1.10 and run `Pkg.resolve()`.
- Prefer optional extensions for units, coordinate objects, file containers, and specialized observatory functionality.

## Testing workflow

- Read the relevant implementation and existing tests before editing.
- Run the narrowest relevant test file or testset first.
- Add regression tests for concrete failure modes, boundary conventions, dimensions, inverses, dependencies, inference, and allocations.
- Warm up code before checking `@allocated`.
- Use a typed helper around `@allocated` when testing Julia 1.10; global or freshly defined testset wrappers can introduce measurement artifacts.
- Cross-check analytical inverses with forward round trips.
- Test zero-dimensional and matrix forms when the API supports them.
- For metadata features, test programmatic construction and tagged parsing separately.
- Add Python-produced or Python-validated fixtures for Python GWCS compatibility features.
- Run `Pkg.test()` only when the change is coherent and targeted tests pass.
- Before completion, run the full suite on current stable Julia and Julia 1.10.
- Run `git diff --check`.
- Use a long timeout when running full package tests (e.g., 10 min).

Prefer the persistent Julia MCP session with the repository environment so Revise can amortize compilation. Do not repeatedly launch fresh Julia processes for normal edit-test loops.

Python GWCS reference work uses:

- Python: `~/.venvs/311_main/bin/python`
- GWCS source: `~/Development/python/gwcs`
- Astropy source: `~/Development/python/astropy`

## Benchmarks and validation

- Add benchmark coverage for new transform families or dispatch mechanisms.
- Compare equivalent FITSWCS.jl and GWCS.jl paths where possible.
- Report allocations as well as time.
- Keep benchmark construction separate from evaluation benchmarks.
- Extend validation scripts when implementing externally defined semantics.
- Document deliberate differences from Python GWCS rather than allowing accidental divergence.

## Documentation and feature tracking

- Update the README for new public APIs, conventions, restrictions, and compatibility differences.
- Update `feature_list.md` when completing or splitting checklist items.
- Keep deferred work explicit, especially when correctness depends on future geometry, units, batching, or serialization design.
