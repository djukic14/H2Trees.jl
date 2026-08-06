module TestSimpleHybridTree

using Test

@testset verbose = true "Simple Hybrid Tree" begin
    include("trees/test_simplehybridtree.jl")
    include("plans/test_splitting.jl")
end

end # module TestSimpleHybridTree
