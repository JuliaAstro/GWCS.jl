# Development

```@meta
CurrentModule = GWCS
```

Run the package tests, which include the doctests, with:

```sh
julia --project -e 'using Pkg; Pkg.test()'
```

Build this documentation, which needs Julia 1.11 or later, with:

```sh
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## Validation

The scripts in `validation/` compare GWCS.jl with Python GWCS and astropy. They need Julia 1.11 or later and
[uv](https://docs.astral.sh/uv/), which installs the Python packages they use into `validation/.venv`:

```sh
uv sync --project validation
julia --project=validation -e 'using Pkg; Pkg.instantiate()'
```

`GWCS_PYTHON` selects another Python environment with the packages in `validation/pyproject.toml` instead. Each
script runs on its own:

```sh
julia --project=validation validation/validate_asdf_parsers.jl
julia --project=validation validation/compare_asdf_parsers.jl
julia --project=validation validation/compare_arithmetic.jl
julia --project=validation validation/compare_bounding_boxes.jl
julia --project=validation validation/validate_asdf_bounding_boxes.jl
julia --project=validation validation/check_schema_versions.jl
julia --project=validation validation/compare_written_asdf.jl
julia --project=validation validation/compare_rewrites.jl
julia --project=validation validation/compare_fits_asdf.jl
julia --project=validation validation/audit_inference.jl [FILE.asdf...]
julia --project=validation validation/compare_roman.jl FILE...
```

`compare_roman.jl` compares GWCS.jl with Python GWCS on Roman WFI level 2 (calibrated image) ASDF files, given as
arguments, whose WCS is the `roman.meta.wcs` node. Python reads them with the `asdf` and `gwcs` packages already in the
validation environment; Roman's own packages, such as `roman_datamodels`, are not needed.

- The parser validation checks ASDF.jl lazy-array orientation and compares identity, affine, 2D rotation, TAN,
  native/celestial rotation, Euler rotation, spherical rotation-sequence, arithmetic compound, bounding-box,
  label-mapper, and region-selector behavior against astropy and Python GWCS.
- `audit_inference.jl` checks that scalar `pixel_to_world`, `world_to_pixel`, `pixel_to_world_jacobian`,
  `world_to_pixel_jacobian`, `in_domain`, and `in_image` calls are inferred and allocation-free, and that matrix calls
  do not allocate per column, across every FITS projection, the supported distortions, spectral, time, tabular, and
  mixed-axis headers, the ASDF fixtures, and programmatic pipelines. ASDF files passed on the command line, such as
  Roman L2 products, are audited too, with and without their explicit inverses.
- `compare_written_asdf.jl` writes the canonical forms that `compose` and `parallel` build, loads them with Python GWCS,
  and compares the results.
- `compare_rewrites.jl` writes back the Python-written quantity, projection, frame axis-order, and metadata-fidelity
  fixtures and checks that Python reads them with the same units, frames, parameter constraints, and results.
- `compare_fits_asdf.jl` writes FITS-derived WCSs, including SIP, TPV, and TPD distortions, to ASDF and checks that
  Python GWCS reads the same coordinates as GWCS.jl and the celestial frames, reference times, and coordinates astropy
  derives from the same headers.
- `check_schema_versions.jl` verifies that the built-in exact-version table matches the installed ASDF schema
  packages.
- The Roman validation compares four detector positions per file in both directions against Python GWCS, within
  `1e-10` degrees forward and `1e-8` pixels inverse.
