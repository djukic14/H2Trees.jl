# `@inferred` checks: concrete-type inference for construction and iterator consumption.
#
# `@inferred` catches a return type that resolves to an abstract/`Union` type, which is the class
# of issue most directly responsible for allocations and dynamic dispatch on a hot path.

using Test
using H2Trees
using StaticArrays
using Graphs

@testset "Tree construction inference" begin
    for N in (1, 2, 3)
        pts = perf_points(N)
        tree = @inferred perf_buildtwontree(pts)
        @test tree isa H2Trees.TwoNTree
    end

    pts = perf_points(3)
    trial = perf_trial_points(3)
    block = @inferred perf_buildblocktree(pts, trial)
    @test block isa H2Trees.BlockTree

    hybrid = @inferred perf_buildsimplehybridtree(pts)
    @test hybrid isa H2Trees.SimpleHybridTree

    # `BoundingBallTreeBuilder.balldata` (default `BoundingBallData`) is a runtime function value
    # carried in the builder struct, not encoded in a type parameter, so `_splitboundingballnode!`
    # cannot fully resolve the constructed node data type at compile time. This is a genuine,
    # pre-existing dispatch gap shared by every `BoundingBallTree`-based construction
    # (`BoundingBallTree`/`KMeansTree`/`MetisTree`/`MetisForest`)
    # and left as a tracked follow-up rather than silently ignored.
    g = perf_graph(length(pts))
    w = perf_weights(length(pts))
    gf = perf_forest_graph(100, 100)
    wf = perf_weights(200)
    @test_broken (@inferred perf_buildboundingballtree(pts)) isa H2Trees.BoundingBallTree
    @test_broken (@inferred perf_buildkmeanstree(pts)) isa H2Trees.BoundingBallTree
    @test_broken (@inferred perf_buildmetistree(pts, g, w)) isa H2Trees.BoundingBallTree
    @test_broken (@inferred perf_buildmetisforest(pts, gf, wf)) isa H2Trees.Forest
end

@testset "Iterator inference" begin
    tree = perf_buildtwontree(perf_points(3))
    rootnode = H2Trees.root(tree)
    leaf = first(H2Trees.leaves(tree))
    lvl = H2Trees.levels(tree)[1]

    @test (@inferred H2Trees.DepthFirstIterator(tree)) isa H2Trees.DepthFirstIterator
    @test (@inferred H2Trees.ChildIterator(tree, rootnode)) isa H2Trees.ChildIterator
    @test (@inferred H2Trees.ParentUpwardsIterator(tree, leaf)) isa
        H2Trees.ParentUpwardsIterator
    @inferred H2Trees.LevelIterator(tree, lvl)
    @inferred H2Trees.NearNodeIterator(tree, rootnode)
    @inferred H2Trees.FarNodeIterator(tree, rootnode)

    @test (@inferred perf_depthfirstsum(tree)) isa Int
    @test (@inferred perf_childrensum(tree)) isa Int
    @test (@inferred perf_parentupwardssum(tree)) isa Int
    @test (@inferred perf_levelsum(tree)) isa Int
    @test (@inferred perf_nearfarsum(tree)) isa Tuple{Int,Int}
end

@testset "Plan / checkadmissibility inference" begin
    plantree = perf_buildplantree(perf_points(3))
    galerkinplans = @inferred perf_buildgalerkinplans(plantree)
    @test galerkinplans isa H2Trees.PlanSet
    report = @inferred perf_checkadmissibility(plantree, galerkinplans)
    @test report isa H2Trees.AdmissibilityReport

    block = perf_buildplanblocktree(perf_points(3), perf_trial_points(3))
    petrovplans = @inferred perf_buildpetrovplans(block)
    @test petrovplans isa H2Trees.PlanSet
    reportpetrov = @inferred perf_checkadmissibility(block, petrovplans)
    @test reportpetrov isa H2Trees.AdmissibilityReport
end
