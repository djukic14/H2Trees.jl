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

@testset "builder validation rejects inconsistent settings" begin
    points = [
        SVector(0.0, 0.0, 0.0),
        SVector(1.0, 1.0, 1.0),
        SVector(0.25, 0.25, 0.25),
        SVector(0.75, 0.75, 0.75),
    ]

    # A `BlockTree`'s two sides share one level scale, so mismatched `minhalfsize` or `root`
    # would silently produce sides that cannot be compared level for level. Both guards sit on
    # single `||` lines, so line coverage from the passing case hides them entirely.
    @test_throws ArgumentError BlockTreeBuilder(;
        test=TwoNTreeBuilder(; minhalfsize=0.25), trial=TwoNTreeBuilder(; minhalfsize=0.5)
    )
    @test_throws ArgumentError BlockTreeBuilder(;
        test=TwoNTreeBuilder(; root=1), trial=TwoNTreeBuilder(; root=2)
    )
    # Matching settings are accepted, so the guards are not rejecting everything.
    @test BlockTreeBuilder(;
        test=TwoNTreeBuilder(; minhalfsize=0.25, root=3),
        trial=TwoNTreeBuilder(; minhalfsize=0.25, root=3),
    ) isa BlockTreeBuilder

    # `minlevel` accepts an `Int` or `AutoMinLevel()`; anything else must be rejected where it
    # is resolved rather than propagating as a confusing failure deeper in construction.
    @test H2Trees._resolve_minlevel(AutoMinLevel(), 4) == 4
    @test H2Trees._resolve_minlevel(2, 4) == 2
    @test_throws ArgumentError H2Trees._resolve_minlevel(2.5, 4)
    @test_throws ArgumentError H2Trees._resolve_minlevel(nothing, 4)

    # A ball-tree splitter is called through `_callsplitwrapper`, which accepts either the
    # `(points, values, level, numsplits)` or the older `(points, values, numsplits)` shape.
    # One that matches neither must say so, instead of surfacing as a `MethodError` from
    # somewhere inside the recursion.
    wrongarity(points, globalpointids) = ([globalpointids], [SVector(0.0, 0.0, 0.0)], [1.0])
    @test_throws ArgumentError buildtree(
        points;
        builder=BoundingBallTreeBuilder(; splitter=wrongarity, numsplits=2, minvalues=1),
    )

    # The older three-argument splitter shape is an explicit BoundingBallTree compatibility
    # allowance for existing custom splitters; the check above is only about rejecting an
    # unsupported arity.
    oldshape(points, globalpointids, numsplits) = (
        [globalpointids[1:2], globalpointids[3:end]],
        [SVector(0.0, 0.0, 0.0), SVector(1.0, 1.0, 1.0)],
        [0.5, 0.5],
    )
    @test buildtree(
        points;
        builder=BoundingBallTreeBuilder(; splitter=oldshape, numsplits=2, minvalues=1),
    ) isa H2Trees.BoundingBallTree
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
