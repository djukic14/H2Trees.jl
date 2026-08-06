using Test
using StaticArrays
using CompScienceMeshes, BEAST
using H2Trees

# `BlockTree`'s `testtree`/`trialtree` are each just a `TwoNTree` going through the same
# `bulkbuildtwontree` pipeline as a standalone `TwoNTree` (see `BlockTree.jl`'s `_blocktwontree`),
# but with root geometry/minlevel resolved JOINTLY across both sides by
# `adjusttwontreeblocktreeparameters` rather than independently. This checks that resolution feeds
# `bulkbuildtwontree` correctly, across the shapes Phase 7 calls out: symmetric Galerkin-like
# (identical point sets), separated Petrov (distinct point sets), asymmetric depth (root sizes
# differ enough that one side's minlevel starts offset into the other's levels), and a real
# BEAST-space geometry pair.

function _detpoints(N::Int, n::Int; offset=0.0, scale=1.0)
    golden = (sqrt(5.0) - 1.0) / 2.0
    return [
        SVector(ntuple(d -> offset + scale * mod(i * golden^d, 1.0) * 10.0, N)) for i in 1:n
    ]
end

function _isomorphic(tree1, node1, tree2, node2)
    d1, d2 = H2Trees.data(tree1, node1), H2Trees.data(tree2, node2)
    (
        d1.sector == d2.sector &&
        d1.center == d2.center &&
        d1.halfsize == d2.halfsize &&
        d1.level == d2.level
    ) || return false
    H2Trees.isleaf(tree1, node1) == H2Trees.isleaf(tree2, node2) || return false
    H2Trees.isleaf(tree1, node1) && return d1.values == d2.values

    children1 = collect(H2Trees.children(tree1, node1))
    children2 = collect(H2Trees.children(tree2, node2))
    length(children1) == length(children2) || return false
    for (c1, c2) in zip(children1, children2)
        _isomorphic(tree1, c1, tree2, c2) || return false
    end
    return true
end

function _samenodes(tree1, tree2)
    return _isomorphic(tree1, H2Trees.root(tree1), tree2, H2Trees.root(tree2))
end

# Rebuilds what each side of a `BlockTree(...)` call should produce, by independently
# re-deriving the resolved root geometry/minlevel via the same public
# `adjusttwontreeblocktreeparameters` `_blocktwontree` calls, and feeding those straight into
# `bulkbuildtwontree` directly (bypassing `_blocktwontree`/`TwoNTree(...)`'s own wiring entirely).
# This checks `_blocktwontree`'s resolution/wiring specifically -- if it ever passed the wrong
# center/halfsize/minlevel/root to one side, this independently-computed reference would diverge
# from what `BlockTree(...)` actually built, even though both ultimately call
# `bulkbuildtwontree`.
function _referenceblocksides(testpositions, trialpositions, builder)
    testcenter, testhalfsize = H2Trees.boundingbox(testpositions)
    trialcenter, trialhalfsize = H2Trees.boundingbox(trialpositions)
    minhalfsize = oftype(testhalfsize, builder.test.minhalfsize)
    _, testroothalfsize, testminlevel, trialroothalfsize, trialminlevel = H2Trees.adjusttwontreeblocktreeparameters(
        testhalfsize, trialhalfsize, minhalfsize
    )

    reftesttree = H2Trees.bulkbuildtwontree(
        testpositions,
        SVector(testcenter...),
        testroothalfsize,
        minhalfsize,
        testminlevel,
        builder.test.root,
        builder.test.minvalues,
        builder.test.protrusion,
    )
    reftrialtree = H2Trees.bulkbuildtwontree(
        trialpositions,
        SVector(trialcenter...),
        trialroothalfsize,
        minhalfsize,
        trialminlevel,
        builder.trial.root,
        builder.trial.minvalues,
        builder.trial.protrusion,
    )
    return reftesttree, reftrialtree, testminlevel, trialminlevel
end

@testset "BlockTree bulk-build equivalence" begin
    @testset "symmetric Galerkin-like: identical test/trial point sets" begin
        points = _detpoints(3, 90)
        builder = H2Trees.BlockTreeBuilder(;
            test=H2Trees.TwoNTreeBuilder(;
                minvalues=5, protrusion=H2Trees.NoProtrusionCheck()
            ),
            trial=H2Trees.TwoNTreeBuilder(;
                minvalues=5, protrusion=H2Trees.NoProtrusionCheck()
            ),
        )
        reftest, reftrial, testminlevel, trialminlevel = _referenceblocksides(
            points, points, builder
        )
        @test testminlevel == trialminlevel == 1

        tree = H2Trees.BlockTree(points, points; builder=builder)
        @test _samenodes(H2Trees.testtree(tree), reftest)
        @test _samenodes(H2Trees.trialtree(tree), reftrial)
    end

    @testset "separated Petrov: distinct, differently-shaped point sets" begin
        testpoints = _detpoints(3, 90)
        trialpoints = _detpoints(3, 60; offset=3.0, scale=0.7)
        builder = H2Trees.BlockTreeBuilder(;
            test=H2Trees.TwoNTreeBuilder(;
                minvalues=7, protrusion=H2Trees.NoProtrusionCheck()
            ),
            trial=H2Trees.TwoNTreeBuilder(;
                minvalues=3, protrusion=H2Trees.NoProtrusionCheck()
            ),
        )
        reftest, reftrial, = _referenceblocksides(testpoints, trialpoints, builder)

        tree = H2Trees.BlockTree(testpoints, trialpoints; builder=builder)
        @test _samenodes(H2Trees.testtree(tree), reftest)
        @test _samenodes(H2Trees.trialtree(tree), reftrial)
    end

    @testset "asymmetric depth: root sizes differ enough that one side starts offset" begin
        testpoints = _detpoints(3, 90; scale=8.0)
        trialpoints = _detpoints(3, 70; offset=1.0, scale=1.0)
        builder = H2Trees.BlockTreeBuilder(;
            test=H2Trees.TwoNTreeBuilder(;
                minvalues=4, protrusion=H2Trees.NoProtrusionCheck()
            ),
            trial=H2Trees.TwoNTreeBuilder(;
                minvalues=4, protrusion=H2Trees.NoProtrusionCheck()
            ),
        )
        reftest, reftrial, testminlevel, trialminlevel = _referenceblocksides(
            testpoints, trialpoints, builder
        )
        # otherwise this case isn't actually exercising asymmetric depth at all
        @test testminlevel != trialminlevel

        tree = H2Trees.BlockTree(testpoints, trialpoints; builder=builder)
        @test _samenodes(H2Trees.testtree(tree), reftest)
        @test _samenodes(H2Trees.trialtree(tree), reftrial)
    end

    @testset "BEAST-space constructor path, real geometry" begin
        mx = CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere5.in")
        )
        my = CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter7.in")
        )
        X = raviartthomas(mx)
        Y = raviartthomas(my)
        testpoints = BEAST.positions(X)
        trialpoints = BEAST.positions(Y)

        builder = H2Trees.BlockTreeBuilder(;
            test=H2Trees.TwoNTreeBuilder(;
                minhalfsize=0.1,
                minvalues=10,
                protrusion=H2Trees.ProtrusionCheck(;
                    max=0.3, compute=H2Trees.BEASTProtrusionFunctor(X)
                ),
            ),
            trial=H2Trees.TwoNTreeBuilder(;
                minhalfsize=0.1,
                minvalues=3,
                protrusion=H2Trees.ProtrusionCheck(;
                    max=0.3, compute=H2Trees.BEASTProtrusionFunctor(Y)
                ),
            ),
        )
        reftest, reftrial, = _referenceblocksides(testpoints, trialpoints, builder)

        tree = H2Trees.BlockTree(testpoints, trialpoints; builder=builder)
        @test _samenodes(H2Trees.testtree(tree), reftest)
        @test _samenodes(H2Trees.trialtree(tree), reftrial)
    end
end
