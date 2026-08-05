module TestPlans

using Test

@testset verbose = true "Plans" begin
    include("test_adjointplans.jl")
    include("test_aggregationplans.jl")
    include("test_disaggregationplans.jl")
    include("test_planapi.jl")
    include("test_checkadmissibility.jl")
end

end # module TestPlans
