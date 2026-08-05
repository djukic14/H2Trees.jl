using H2Trees
using StaticArrays
using Test

@testset "Builder defaults" begin
    twon = TwoNTreeBuilder()
    @test twon.minhalfsize == 0
    @test twon.minvalues == 70
    @test twon.minlevel == AutoMinLevel()
    @test twon.root == 1
    @test twon.protrusion == AutoProtrusionCheck()

    kmeans = KMeansTreeBuilder()
    @test kmeans.numberofclusters == 4
    @test kmeans.minvalues == 70

    metis = MetisTreeBuilder()
    @test metis.numdivisions == 4
    @test metis.minvalues == 4

    forest = MetisForestBuilder()
    @test forest.treebuilder isa MetisTreeBuilder
    @test forest.treebuilder.numdivisions == 4

    block = BlockTreeBuilder()
    @test block.test isa TwoNTreeBuilder
    @test block.trial isa TwoNTreeBuilder
    @test block.test.minhalfsize == block.trial.minhalfsize
end

@testset "buildtree entry points" begin
    points = [
        SVector(0.0, 0.0, 0.0),
        SVector(1.0, 1.0, 1.0),
        SVector(0.25, 0.25, 0.25),
        SVector(0.75, 0.75, 0.75),
    ]

    # Default single-argument form builds a TwoNTree and matches the explicit constructor.
    default = buildtree(points)
    @test default isa H2Trees.TwoNTree
    @test !isempty(H2Trees.leaves(default))

    # The direct constructor also defaults its builder (not just `buildtree`).
    @test TwoNTree(points) isa H2Trees.TwoNTree
    @test H2Trees.numberofnodes(TwoNTree(points)) == H2Trees.numberofnodes(default)

    twonbuilder = TwoNTreeBuilder(; minhalfsize=0.25, minvalues=1)
    viabuild = buildtree(points; builder=twonbuilder)
    viactor = TwoNTree(points; builder=twonbuilder)
    @test H2Trees.numberofnodes(viabuild) == H2Trees.numberofnodes(viactor)
    @test H2Trees.leaves(viabuild) == H2Trees.leaves(viactor)

    # Two point sets build a BlockTree.
    block = buildtree(
        points,
        points;
        builder=BlockTreeBuilder(;
            test=TwoNTreeBuilder(; minhalfsize=0.25, minvalues=1),
            trial=TwoNTreeBuilder(; minhalfsize=0.25, minvalues=2),
        ),
    )
    @test block isa H2Trees.BlockTree
    @test H2Trees.minhalfsize(H2Trees.testtree(block)) ==
        H2Trees.minhalfsize(H2Trees.trialtree(block))

    # The direct BlockTree constructor also defaults its builder, mirroring TwoNTree(points) above.
    @test BlockTree(points, points) isa H2Trees.BlockTree
    blockbuilder = BlockTreeBuilder(;
        test=TwoNTreeBuilder(; minhalfsize=0.25, minvalues=1),
        trial=TwoNTreeBuilder(; minhalfsize=0.25, minvalues=2),
    )
    blockviabuild = buildtree(points, points; builder=blockbuilder)
    blockviactor = BlockTree(points, points; builder=blockbuilder)
    @test H2Trees.numberofnodes(H2Trees.testtree(blockviabuild)) ==
        H2Trees.numberofnodes(H2Trees.testtree(blockviactor))
    @test H2Trees.numberofnodes(H2Trees.trialtree(blockviabuild)) ==
        H2Trees.numberofnodes(H2Trees.trialtree(blockviactor))

    # A BoundingBallTreeBuilder routes through the single-argument entry point too.
    newsplit(points, globalpointids, level, numsplits) = (
        [globalpointids[1:2], globalpointids[3:4]],
        [SVector(level, 0.0, 0.0), SVector(level, 1.0, 1.0)],
        [0.5, 0.5],
    )
    ball = buildtree(
        points;
        builder=BoundingBallTreeBuilder(; splitter=newsplit, numsplits=2, minvalues=1),
    )
    @test ball isa H2Trees.BoundingBallTree

    # An unsupported (points, builder) combination is rejected with a clear error.
    @test_throws ArgumentError buildtree(points; builder=42)
end

@testset "show methods" begin
    b = TwoNTreeBuilder(; minhalfsize=0.25, minvalues=3)
    s = sprint(show, b)
    @test occursin("TwoNTreeBuilder", s)
    @test occursin("minvalues=3", s)
    @test occursin("minlevel=auto", s)
    @test occursin("root=1", s)

    @test occursin("numberofclusters=4", sprint(show, KMeansTreeBuilder()))
    @test occursin("numdivisions=4", sprint(show, MetisTreeBuilder()))
    @test occursin("BlockTreeBuilder", sprint(show, BlockTreeBuilder()))
    @test occursin(
        "hybridhalfsize=0.5", sprint(show, SimpleHybridTreeBuilder(; hybridhalfsize=0.5))
    )

    points = [SVector(Float64(i), 0.0, 0.0) for i in 1:50]
    tree = buildtree(points; builder=TwoNTreeBuilder(; minhalfsize=0.0, minvalues=5))
    ts = sprint(show, tree)
    @test occursin("TwoNTree{3,Float64}", ts)
    @test occursin("nodes=", ts)
    @test occursin("leaves=", ts)
    @test occursin("levels=", ts)
    @test occursin("root=", ts)
end
