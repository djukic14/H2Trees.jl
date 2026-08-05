module TestH2MetisTrees

using Test

@testset verbose = true "H2MetisTrees" begin
    include("test_H2MetisTrees.jl")
    include("test_H2MetisForest.jl")
end

end # module TestH2MetisTrees
