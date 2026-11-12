using Documenter
using Documenter.Remotes: GitHub
using GWCS

# The table of Python equivalents comes from the docstrings' "Python equivalent:" lines, so the two cannot drift apart.
function python_equivalents_page(path)
    rows = String[]
    for name in sort!(filter(!=(:GWCS), names(GWCS)))
        text = string(Docs.doc(getfield(GWCS, name)))
        found = match(r"Python equivalent: ((?:.|\n(?!\n))*)", text)
        found === nothing && continue
        # A pipe inside a cell, as in the `|` operator, would end it.
        push!(rows, "| [`$name`](@ref) | $(replace(found.captures[1], '\n' => ' ', '|' => "\\|")) |")
    end
    header = """
    # Python equivalents

    ```@meta
    CurrentModule = GWCS
    EditURL = "../make.jl"
    ```

    Each GWCS.jl name with a counterpart in Python GWCS or astropy, as its docstring states it.

    | GWCS.jl | Python |
    |:--------|:-------|
    """
    write(path, header * join(rows, "\n") * "\n")
end

python_equivalents_page(joinpath(@__DIR__, "src", "python_equivalents.md"))

# Pretty URLs need a web server; local builds use plain .html links that open from disk.
const CI = get(ENV, "CI", "false") == "true"

makedocs(;
    modules=[GWCS],
    authors="Chris Garling",
    repo=GitHub("cgarling/GWCS.jl"),
    sitename="GWCS.jl",
    format=Documenter.HTML(; prettyurls=CI, canonical="https://cgarling.github.io/GWCS.jl/stable/", edit_link="main"),
    pages=[
        "Home" => "index.md",
        "Building and evaluating a WCS" => "pipelines.md",
        "Transforms" => "transforms.md",
        "FITS WCS" => "fits.md",
        "ASDF and metadata" => "asdf.md",
        "Conventions and differences" => "conventions.md",
        "Python equivalents" => "python_equivalents.md",
        "API reference" => [
            "WCS" => "api/wcs.md",
            "Frames and domains" => "api/frames.md",
            "Transforms" => "api/transforms.md",
            "Metadata and extensions" => "api/metadata.md",
        ],
        "Development" => "development.md",
    ],
    # Doctests run with the test suite.
    doctest=false,
    checkdocs=:public,
)

deploydocs(; repo="github.com/cgarling/GWCS.jl", devbranch="main", push_preview=true, versions=["stable" => "v^", "v#.#", "dev" => "dev"])
