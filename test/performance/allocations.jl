# Allocation-budget and zero-allocation hot-path checks. See `budgets.jl` for the documented
# per-family constants and the measured baselines they were set from.

using Test
using H2Trees
using StaticArrays
using Graphs

# Warms up (compiles + runs twice, discarding results) before measuring so `@timed`/`@allocated`
# never counts JIT compilation, which otherwise dwarfs the real allocation and makes every budget
# meaningless.
function perf_measure(f, args...)
    f(args...)
    f(args...)
    stats = @timed f(args...)
    return stats.bytes, stats.value
end

function perf_allocationratio(f, args...)
    bytes, value = perf_measure(f, args...)
    return bytes / Base.summarysize(value)
end

@testset "Construction allocation budgets" begin
    for N in (1, 2, 3)
        for pts in (perf_points(N), perf_points(N; large=true))
            @test perf_allocationratio(perf_buildtwontree, pts) <=
                ALLOCATION_RATIO_BUDGET[:TwoNTree]
        end
    end

    for (pts, trial) in (
        (perf_points(3), perf_trial_points(3)),
        (perf_points(3; large=true), perf_trial_points(3; large=true)),
    )
        @test perf_allocationratio(perf_buildblocktree, pts, trial) <=
            ALLOCATION_RATIO_BUDGET[:BlockTree]
    end

    for pts in (perf_points(3), perf_points(3; large=true))
        @test perf_allocationratio(perf_buildsimplehybridtree, pts) <=
            ALLOCATION_RATIO_BUDGET[:SimpleHybridTree]
        @test perf_allocationratio(perf_buildboundingballtree, pts) <=
            ALLOCATION_RATIO_BUDGET[:BoundingBallTree]
        @test perf_allocationratio(perf_buildkmeanstree, pts) <=
            ALLOCATION_RATIO_BUDGET[:KMeansTree]
    end

    for pts in (perf_points(3), perf_points(3; large=true))
        g = perf_graph(length(pts))
        w = perf_weights(length(pts))
        @test perf_allocationratio(perf_buildmetistree, pts, g, w) <=
            ALLOCATION_RATIO_BUDGET[:MetisTree]
    end

    let pts = perf_points(3), gf = perf_forest_graph(100, 100), wf = perf_weights(200)
        @test perf_allocationratio(perf_buildmetisforest, pts, gf, wf) <=
            ALLOCATION_RATIO_BUDGET[:MetisForest]
    end
end

@testset "Construction allocation scaling" begin
    # A 10x larger point set should not need a much larger allocation-per-final-byte ratio --
    # `ALLOCATION_RATIO_SCALING_TOLERANCE` bounds how much worse `large` is allowed to be than
    # `small`, normalizing away that both grow in absolute bytes as N grows.
    smallpts = perf_points(3)
    largepts = perf_points(3; large=true)

    ratiosmall = perf_allocationratio(perf_buildtwontree, smallpts)
    ratiolarge = perf_allocationratio(perf_buildtwontree, largepts)
    @test ratiolarge <= ratiosmall * ALLOCATION_RATIO_SCALING_TOLERANCE

    ratiosmall = perf_allocationratio(perf_buildboundingballtree, smallpts)
    ratiolarge = perf_allocationratio(perf_buildboundingballtree, largepts)
    @test ratiolarge <= ratiosmall * ALLOCATION_RATIO_SCALING_TOLERANCE

    ratiosmall = perf_allocationratio(perf_buildkmeanstree, smallpts)
    ratiolarge = perf_allocationratio(perf_buildkmeanstree, largepts)
    @test ratiolarge <= ratiosmall * ALLOCATION_RATIO_SCALING_TOLERANCE

    gsmall = perf_graph(length(smallpts))
    wsmall = perf_weights(length(smallpts))
    glarge = perf_graph(length(largepts))
    wlarge = perf_weights(length(largepts))
    ratiosmall = perf_allocationratio(perf_buildmetistree, smallpts, gsmall, wsmall)
    ratiolarge = perf_allocationratio(perf_buildmetistree, largepts, glarge, wlarge)
    @test ratiolarge <= ratiosmall * ALLOCATION_RATIO_SCALING_TOLERANCE
end

@testset "Plan allocation budgets" begin
    plantree = perf_buildplantree(perf_points(3))
    @test perf_allocationratio(perf_buildgalerkinplans, plantree) <=
        ALLOCATION_RATIO_PLAN_BUDGET

    planblock = perf_buildplanblocktree(perf_points(3), perf_trial_points(3))
    @test perf_allocationratio(perf_buildpetrovplans, planblock) <=
        ALLOCATION_RATIO_PLAN_BUDGET
end

@testset "checkadmissibility allocation budget" begin
    # `checkadmissibility`'s natural "result" (`AdmissibilityReport`) is tiny regardless of tree
    # size, so unlike construction/plan budgets above, the ratio here is measured against
    # `summarysize(tree)` -- see `budgets.jl` for why.
    plantree = perf_buildplantree(perf_points(3))
    galerkinplans = perf_buildgalerkinplans(plantree)
    bytes, _ = perf_measure(perf_checkadmissibility, plantree, galerkinplans)
    @test bytes / Base.summarysize(plantree) <= ALLOCATION_RATIO_CHECKADMISSIBILITY_BUDGET

    planblock = perf_buildplanblocktree(perf_points(3), perf_trial_points(3))
    petrovplans = perf_buildpetrovplans(planblock)
    bytes, _ = perf_measure(perf_checkadmissibility, planblock, petrovplans)
    @test bytes / Base.summarysize(planblock) <= ALLOCATION_RATIO_CHECKADMISSIBILITY_BUDGET
end

@testset "Zero-allocation iterator hot paths" begin
    tree = perf_buildtwontree(perf_points(3))
    rootnode = H2Trees.root(tree)
    leaf = first(H2Trees.leaves(tree))
    lvl = H2Trees.levels(tree)[end]

    # `ChildIterator`/`ParentUpwardsIterator`/`LevelIterator` consume existing tree structure
    # (linked-list-style sibling pointers, or a cached level vector) without allocating anything
    # new -- these are the iterators the plan's "zero allocations" bar actually applies to.
    function sumchildren(tree, node)
        s = 0
        for child in H2Trees.ChildIterator(tree, node)
            s += child
        end
        return s
    end
    sumchildren(tree, rootnode)
    @test (@allocated sumchildren(tree, rootnode)) == 0

    function sumparents(tree, node)
        s = 0
        for p in H2Trees.ParentUpwardsIterator(tree, node)
            s += p
        end
        return s
    end
    sumparents(tree, leaf)
    @test (@allocated sumparents(tree, leaf)) == 0

    function sumlevel(tree, level)
        s = 0
        for n in H2Trees.LevelIterator(tree, level)
            s += n
        end
        return s
    end
    sumlevel(tree, lvl)
    @test (@allocated sumlevel(tree, lvl)) == 0

    # `numberofvalues` recurses over `children(tree, node)` (like `appendvalues!`/`foreachvalue`/
    # `anyvalue`) rather than materializing a leaf list via `leaves(tree, node)` -- covers the
    # root, an internal non-root node, and a leaf, since the root case used to be cheap already
    # (`leaves(tree, root(tree))` is a cached-vector `copy`) while a non-root internal node used to
    # pay for a fresh `DepthFirstIterator` + `collect`.
    plantree = perf_buildplantree(perf_points(3))
    internalnode = 0
    for node in H2Trees.DepthFirstIterator(plantree)
        if node != H2Trees.root(plantree) && !H2Trees.isleaf(plantree, node)
            internalnode = node
            break
        end
    end
    @test internalnode != 0
    planleaf = first(H2Trees.leaves(plantree))

    H2Trees.numberofvalues(plantree, H2Trees.root(plantree))
    @test (@allocated H2Trees.numberofvalues(plantree, H2Trees.root(plantree))) == 0
    H2Trees.numberofvalues(plantree, internalnode)
    @test (@allocated H2Trees.numberofvalues(plantree, internalnode)) == 0
    H2Trees.numberofvalues(plantree, planleaf)
    @test (@allocated H2Trees.numberofvalues(plantree, planleaf)) == 0
end

@testset "Iterator consumption allocations stay small" begin
    # `DepthFirstIterator`/`NearNodeIterator`/`FarNodeIterator` are NOT zero-allocation by design
    # -- `DepthFirstIterator` builds a traversal stack, and the near/far iterators collect a
    # filtered node list -- each fresh call legitimately allocates that one-time setup. This just
    # bounds it to a small constant per node visited, so it can't silently regress into something
    # that scales badly with tree size.
    tree = perf_buildtwontree(perf_points(3))
    # A small tree's traversal-stack/filtered-list setup cost dominates per-node bytes (measured
    # ~154 bytes/node at 5 nodes, ~65 bytes/node at 67 nodes) -- the `max(nnodes, 10)` floor keeps
    # the budget meaningful for a tree this shallow instead of shrinking below the fixed setup
    # cost.
    nnodes = max(H2Trees.numberofnodes(tree), 10)

    bytes, _ = perf_measure(perf_depthfirstsum, tree)
    @test bytes <= 200 * nnodes

    bytes, _ = perf_measure(perf_childrensum, tree)
    @test bytes <= 200 * nnodes

    bytes, _ = perf_measure(perf_nearfarsum, tree)
    @test bytes <= 200 * nnodes
end
