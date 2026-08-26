using ForwardDiff
using GWCS
using StaticArrays
using Test

@testset "GWCS" begin
    include("test_transforms.jl")
    include("test_callable.jl")
    include("test_composition.jl")
    include("test_pipeline.jl")
    include("test_api.jl")
    include("test_numerical_inverse.jl")
    include("test_metadata.jl")
    include("test_registration.jl")
end
