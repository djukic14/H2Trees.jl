module TestSEBB

using Test
using H2Trees
using StaticArrays
using LinearAlgebra
using Random
using BoundingSphere

@testset verbose = true "SEBB" begin
    include("test_helpers.jl")
    include("test_ball.jl")
    include("test_roots.jl")
    include("test_known_solutions.jl")
    include("test_constructed_optima.jl")
    include("test_certificates.jl")
    include("test_metamorphic.jl")
    include("test_degeneracies.jl")
    include("test_fallback.jl")
    include("test_high_precision.jl")
    include("test_randomized_1d.jl")
    include("test_inference.jl")
    include("test_allocations.jl")
end

end # module TestSEBB
