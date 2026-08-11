module TestBlockBallTree

# A `BlockTree` whose two sides are `BoundingBallTree`s.
#
# `BlockTree(testpositions, trialpositions; builder)` only ever builds `TwoNTree` sides, but the
# two-argument `BlockTree(testcluster, trialcluster)` accepts any pair of trees, and the package
# carries trait methods written specifically for this shape: the two-tree
# `isnear(..., ::isBoundingBallTree, ::isBoundingBallTree)` predicate and
# `nodegap(..., ::isBoundingBallTree)`. Nothing in the suite built that combination, so those
# methods were only ever reached through the single-tree path: which is how a `kwargs`
# splatted positionally instead of by keyword survived in the two-tree one: with no keywords the
# splat is empty and the call happens to work.

using H2Trees
using StaticArrays
using LinearAlgebra
using BoundingSphere
using Random
using Test

# Binary median splitter (same shape as the one in `test_boundingballtree_sebb.jl`).
function median_split(points, ids, level, numsplits)
    dims = length(points[1])
    coords = [points[i] for i in ids]
    spreads = [
        maximum(p[d] for p in coords) - minimum(p[d] for p in coords) for d in 1:dims
    ]
    d = argmax(spreads)
    sorted = ids[sortperm([points[i][d] for i in ids])]
    mid = cld(length(ids), 2)
    parts = [sorted[1:mid], sorted[(mid + 1):end]]
    centers = SVector{dims,Float64}[]
    radii = Float64[]
    for part in parts
        c, r = BoundingSphere.boundingsphere([points[i] for i in part])
        push!(centers, SVector{dims}(c))
        push!(radii, r)
    end
    return parts, centers, radii
end

function ballblocktree(seed=1; separation=4.0, n=80)
    rng = MersenneTwister(seed)
    testpoints = [SVector(randn(rng), randn(rng), randn(rng)) for _ in 1:n]
    trialpoints = [SVector(randn(rng) + separation, randn(rng), randn(rng)) for _ in 1:n]
    builder = BoundingBallTreeBuilder(; splitter=median_split, numsplits=2, minvalues=8)
    return H2Trees.BlockTree(
        buildtree(testpoints; builder=builder), buildtree(trialpoints; builder=builder)
    )
end

# Independent restatement of `isnearradius`, so the expectation is not the implementation.
function refnear(testtree, trialtree, testnode, trialnode, η)
    ca = H2Trees.center(testtree, testnode)
    cb = H2Trees.center(trialtree, trialnode)
    ra = H2Trees.radius(testtree, testnode)
    rb = H2Trees.radius(trialtree, trialnode)
    d = norm(ca - cb)
    (d + rb <= ra || d + ra <= rb) && return true
    return d <= η * (1 + 100 * eps(Float64)) * (ra + rb)
end

@testset "two-tree ball predicate honours its keywords (regression)" begin
    tree = ballblocktree()
    testtree = H2Trees.testtree(tree)
    trialtree = H2Trees.trialtree(tree)
    @test H2Trees.treetrait(tree) isa H2Trees.isBlockTree
    @test H2Trees.treetrait(testtree) isa H2Trees.isBoundingBallTree

    # The keyword travelled into the predicate as a POSITIONAL argument, so every call carrying
    # one raised `MethodError: no method matching isnearradius(..., ::Pair{Symbol,Float64})`.
    # The default call kept working, which is why this went unnoticed.
    for η in (1.0, 1.5, 4.0)
        predicate = H2Trees.isnear(; η=η)(tree)
        for testnode in H2Trees.leaves(testtree), trialnode in H2Trees.leaves(trialtree)
            @test predicate(testtree, trialtree, testnode, trialnode) ==
                refnear(testtree, trialtree, testnode, trialnode, η)
        end
    end
end

@testset "a larger eta can only enlarge the near set" begin
    tree = ballblocktree()
    testtree = H2Trees.testtree(tree)
    trialtree = H2Trees.trialtree(tree)
    pairs = [(a, b) for a in H2Trees.leaves(testtree) for b in H2Trees.leaves(trialtree)]
    counts = Int[]
    for η in (1.0, 2.0, 8.0)
        predicate = H2Trees.isnear(; η=η)(tree)
        push!(counts, count(p -> predicate(testtree, trialtree, p...), pairs))
    end
    @test issorted(counts)
    # Not a vacuous check: the margin has to actually change something between the extremes,
    # otherwise the equality above would hold for a predicate that ignored `η` entirely.
    @test first(counts) < last(counts)
end

@testset "eta below one is rejected on the two-tree path too" begin
    tree = ballblocktree()
    testtree = H2Trees.testtree(tree)
    trialtree = H2Trees.trialtree(tree)
    predicate = H2Trees.isnear(; η=0.5)(tree)
    @test_throws ArgumentError predicate(
        testtree, trialtree, H2Trees.root(testtree), H2Trees.root(trialtree)
    )
end

@testset "near/far plans on a ball block tree" begin
    # The configuration works end to end; this pins that, and pins that a widened margin
    # reaches `checkadmissibility` rather than erroring on the way in.
    tree = ballblocktree(2)
    @test !isempty(H2Trees.nearinteractions(tree))

    plans = buildplans(tree)
    @test H2Trees.ispetrov(plans)
    @test H2Trees.checkadmissibility(tree, plans; throw=false).ok

    widened = H2Trees.checkadmissibility(
        tree, plans; isnear=H2Trees.isnear(; η=1.0), throw=false
    )
    @test widened isa H2Trees.AdmissibilityReport
end

end # module TestBlockBallTree
