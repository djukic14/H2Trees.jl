using Test, TestItems, TestItemRunner

@testitem "Iterators" begin
    include("test_iterators.jl")
end

@testitem "Near interactions" begin
    include("test_nearinteractions.jl")
end

@testitem "TwoNTree" begin
    include("test_TwoNTree.jl")
end

@testitem "Plans" begin
    include("plans/test_adjointplans.jl")
    include("plans/test_aggregationplans.jl")
    include("plans/test_disaggregationplans.jl")
    include("plans/test_planapi.jl")
end

@testitem "H2NURBSTrees" begin
    using NURBS
    using BEAST
    BEASTnurbs = Base.get_extension(BEAST, :BEASTnurbs)
end

@testitem "H2BEASTTrees" begin
    include("H2BEASTTrees/test_quadpointstree.jl")
end

@testitem "Translations" begin
    include("translations/test_translations.jl")
end

@testitem "Plots" begin
    include("H2PlotlyJSTrees/test_plots.jl")
end

# @testitem "Code quality (Aqua.jl)" begin
#     using Aqua
#     using H2Trees
#     Aqua.test_all(H2Trees)
# end
# @testitem "Code linting (JET.jl)" begin
#     using JET
#     using H2Trees
#     JET.test_package(H2Trees; target_defined_modules=true)
# end

@testitem "Code formatting (JuliaFormatter.jl)" begin
    using JuliaFormatter
    using H2Trees
    pkgpath = pkgdir(H2Trees)
    @test JuliaFormatter.format(pkgpath, overwrite=false)
end

@run_package_tests verbose = true
