module TestH2BEASTTrees

using Test

@testset verbose = true "H2BEASTTrees" begin
    include("test_operators.jl")
    include("test_adjacencygraphs.jl")
end

end # module TestH2BEASTTrees
