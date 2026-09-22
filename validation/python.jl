function validation_python()
    local_python = joinpath(@__DIR__, ".venv", Sys.iswindows() ? "Scripts" : "bin", Sys.iswindows() ? "python.exe" : "python")
    python = get(ENV, "GWCS_PYTHON", isfile(local_python) ? local_python : get(ENV, "PYTHON", local_python))
    isfile(python) || error("validation Python environment not found; run `uv sync --project validation` or set GWCS_PYTHON or PYTHON")
    return python
end
