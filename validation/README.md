# Validation

Create the pinned Python validation environment from the repository root:

```sh
uv sync --project validation
```

Julia validation scripts use `validation/.venv/bin/python` by default. Set `GWCS_PYTHON` to override it.

Run a validation script with its Julia environment, for example:

```sh
julia --project=validation validation/compare_arithmetic.jl
julia --project=validation validation/check_schema_versions.jl
```
