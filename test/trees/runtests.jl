module TestTwoNTree

using Test

@testset verbose = true "TwoNTree" begin
    include("test_TwoNTree.jl")
    include("test_accessors.jl")
    include("test_uniformseparationdepth.jl")
    include("test_checksubdivision.jl")
    include("test_buildtree.jl")
    include("test_bulkbuild.jl")
    include("test_blocktree_bulkbuild.jl")
    include("test_blockballtree.jl")
end

end # module TestTwoNTree
