using Documenter

@testset "Documentation" begin
    # Documenter's checkdocs only checks that existing docstrings appear in the manual, so check that every exported
    # and public name has one.
    if isdefined(Base.Docs, :undocumented_names)  # Julia 1.11 and later
        @test isempty(Base.Docs.undocumented_names(GWCS))
    end

    DocMeta.setdocmeta!(GWCS, :DocTestSetup, :(using GWCS); recursive=true)
    doctest(GWCS; manual=false)
end
