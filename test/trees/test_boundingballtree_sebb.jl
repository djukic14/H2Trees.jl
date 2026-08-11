module TestBoundingBallTreeSEBB

# Phase 6: H2Trees <-> SEBB adapter integration.

using H2Trees
using StaticArrays
using LinearAlgebra
using BoundingSphere
using Random
using Test

const SBT = H2Trees.SEBB

# Recursive median splitter producing a genuine multi-level BoundingBallTree.
function median_split(points, ids, level, numsplits)
    coords = [points[i] for i in ids]
    dims = length(points[1])
    spreads = [
        maximum(p[d] for p in coords) - minimum(p[d] for p in coords) for d in 1:dims
    ]
    d = argmax(spreads)
    order = sortperm([points[i][d] for i in ids])
    sorted = ids[order]
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

# `numsplits`-way splitter, so internal nodes get MORE than two children. The binary
# `median_split` above only ever exercises SEBB's closed-form two-ball case; three or more
# children are what routes through the Gram system and the tangency quadratic.
function nway_split(points, ids, level, numsplits)
    dims = length(points[1])
    sorted = ids[sortperm([points[i][1] for i in ids])]
    n = length(sorted)
    k = min(numsplits, n)
    bounds = [round(Int, n * j / k) for j in 0:k]
    parts = filter(!isempty, [sorted[(bounds[j] + 1):bounds[j + 1]] for j in 1:k])
    centers = SVector{dims,Float64}[]
    radii = Float64[]
    for part in parts
        c, r = BoundingSphere.boundingsphere([points[i] for i in part])
        push!(centers, SVector{dims}(c))
        push!(radii, r)
    end
    return parts, centers, radii
end

# Independent residual (does not use SEBB predicates).
tree_residual(cc, cr, xc, xr) = norm(xc .- cc) + xr - cr

function childballs(tree, node)
    N = length(H2Trees.center(tree, node))
    T = eltype(H2Trees.center(tree, node))
    return [
        SBT.Ball(H2Trees.center(tree, c), H2Trees.radius(tree, c)) for
        c in H2Trees.children(tree, node)
    ]
end

# Test-local sequential pairwise merge (the OLD coarse approximation) for the baseline check.
function pairwise_merge(balls)
    c = SBT.center(balls[1])
    r = SBT.radius(balls[1])
    for k in 2:length(balls)
        c2 = SBT.center(balls[k])
        r2 = SBT.radius(balls[k])
        diff = c .- c2
        dn = norm(diff)
        if dn + r2 <= r
            # keep
        elseif dn + r <= r2
            c, r = c2, r2
        else
            c = 0.5 .* (c .+ c2 .+ (r - r2) .* diff ./ dn)
            r = 0.5 * (r + r2 + dn)
        end
    end
    return c, r
end

@testset "traversal order is child-before-parent (post-order)" begin
    # Regression: locks in the empirically verified fact that DepthFirstIterator on this
    # branch is POST-order. If it ever became pre-order, updateradii! would process parents
    # before children and this test would catch it.
    rng = MersenneTwister(1)
    points = [SVector(randn(rng), randn(rng), randn(rng)) for _ in 1:32]
    tree = buildtree(
        points;
        builder=BoundingBallTreeBuilder(; splitter=median_split, numsplits=2, minvalues=1),
    )
    order = collect(H2Trees.DepthFirstIterator(tree))
    pos = Dict(n => i for (i, n) in enumerate(order))
    for node in order
        p = H2Trees.parent(tree, node)
        if p != 0 && haskey(pos, p)
            @test pos[node] < pos[p]
        end
    end
    @test order[end] == H2Trees.root(tree)   # root last
end

@testset "adapter matches direct SEBB on children" begin
    rng = MersenneTwister(2)
    points = [SVector(randn(rng), randn(rng), randn(rng)) for _ in 1:40]
    tree = buildtree(
        points;
        builder=BoundingBallTreeBuilder(; splitter=median_split, numsplits=2, minvalues=1),
    )
    for node in H2Trees.DepthFirstIterator(tree)
        H2Trees.isleaf(tree, node) && continue
        cbs = childballs(tree, node)
        direct = SBT.smallest_enclosing_ball(cbs)
        ac, ar = H2Trees.boundingsphere(tree, node)
        @test norm(ac .- SBT.center(direct)) <= 1e-9 * max(1.0, ar)
        @test abs(ar - SBT.radius(direct)) <= 1e-9 * max(1.0, ar)
    end
end

@testset "every updated parent contains every child" begin
    rng = MersenneTwister(3)
    points = [SVector(randn(rng), randn(rng), randn(rng)) for _ in 1:50]
    tree = buildtree(
        points;
        builder=BoundingBallTreeBuilder(; splitter=median_split, numsplits=2, minvalues=1),
    )
    for node in H2Trees.DepthFirstIterator(tree)
        H2Trees.isleaf(tree, node) && continue
        cc = H2Trees.center(tree, node)
        cr = H2Trees.radius(tree, node)
        for child in H2Trees.children(tree, node)
            res = tree_residual(
                cc, cr, H2Trees.center(tree, child), H2Trees.radius(tree, child)
            )
            @test res <= 1e-8 * max(1.0, cr)
        end
    end
end

@testset "leaf preservation" begin
    rng = MersenneTwister(4)
    points = [SVector(randn(rng), randn(rng), randn(rng)) for _ in 1:20]
    tree = buildtree(
        points;
        builder=BoundingBallTreeBuilder(; splitter=median_split, numsplits=2, minvalues=1),
    )
    for leaf in H2Trees.leaves(tree)
        c0 = H2Trees.center(tree, leaf)
        r0 = H2Trees.radius(tree, leaf)
        c1, r1 = H2Trees.boundingsphere(tree, leaf)
        @test c1 == c0
        @test r1 == r0
    end
end

@testset "exact radius <= pairwise radius, both enclose children" begin
    rng = MersenneTwister(5)
    points = [SVector(randn(rng), randn(rng), randn(rng)) for _ in 1:60]
    tree = buildtree(
        points;
        builder=BoundingBallTreeBuilder(; splitter=median_split, numsplits=2, minvalues=1),
    )
    for node in H2Trees.DepthFirstIterator(tree)
        H2Trees.isleaf(tree, node) && continue
        cbs = childballs(tree, node)
        exact = SBT.smallest_enclosing_ball(cbs)
        pc, pr = pairwise_merge(cbs)
        @test SBT.radius(exact) <= pr + 1e-9 * max(1.0, pr)
        # both enclose all children
        for b in cbs
            @test tree_residual(
                SBT.center(exact), SBT.radius(exact), SBT.center(b), SBT.radius(b)
            ) <= 1e-8 * max(1.0, SBT.radius(exact))
            @test tree_residual(pc, pr, SBT.center(b), SBT.radius(b)) <= 1e-8 * max(1.0, pr)
        end
    end
end

@testset "bottom-up regression (fails under parent-before-child order)" begin
    # Manual 3-level tree with deliberately-too-small internal/root radii. Only a correct
    # bottom-up update grows node 2 to enclose leaves 4,5 BEFORE the root encloses node 2.
    sv(x, y, z) = SVector{3,Float64}(x, y, z)
    nodes = [
        H2Trees.Node(H2Trees.BoundingBallData(Int[], sv(0, 0, 0), 0.01, 1), 0, 0, 2),  # 1 root
        H2Trees.Node(H2Trees.BoundingBallData(Int[], sv(0, 0, 0), 0.01, 2), 3, 1, 4),  # 2 internal (wrong small r)
        H2Trees.Node(H2Trees.BoundingBallData([3], sv(0, 6, 0), 1.0, 2), 0, 1, 0),      # 3 leaf
        H2Trees.Node(H2Trees.BoundingBallData([1], sv(0, 0, 0), 1.0, 3), 5, 2, 0),      # 4 leaf
        H2Trees.Node(H2Trees.BoundingBallData([2], sv(4, 0, 0), 1.0, 3), 0, 2, 0),      # 5 leaf
    ]
    tree = H2Trees.BoundingBallTree(nodes, 1, sv(0, 0, 0), 0.01, [[1], [2, 3], [4, 5]])
    H2Trees.updateradii!(tree; update=H2Trees.boundingsphere)

    # node 2 must now enclose leaves 4 and 5
    c2 = H2Trees.center(tree, 2)
    r2 = H2Trees.radius(tree, 2)
    for leaf in (4, 5)
        @test tree_residual(
            c2, r2, H2Trees.center(tree, leaf), H2Trees.radius(tree, leaf)
        ) <= 1e-9 * max(1.0, r2)
    end

    # root must enclose all three leaves (only true if node 2 was grown BEFORE the root)
    cr = H2Trees.center(tree, 1)
    rr = H2Trees.radius(tree, 1)
    for leaf in (3, 4, 5)
        @test tree_residual(
            cr, rr, H2Trees.center(tree, leaf), H2Trees.radius(tree, leaf)
        ) <= 1e-9 * max(1.0, rr)
    end
end

@testset "boundingsphereofspheres delegates to exact two-ball SEBB" begin
    c1 = SVector(0.0, 0.0)
    c2 = SVector(4.0, 0.0)
    c, r = H2Trees.boundingsphereofspheres(c1, 1.0, c2, 2.0)
    ref = SBT.smallest_enclosing_ball([SBT.Ball(c1, 1.0), SBT.Ball(c2, 2.0)])
    @test norm(c .- SBT.center(ref)) <= 1e-12
    @test abs(r - SBT.radius(ref)) <= 1e-12
    # nested case returns the containing ball
    cc, rr = H2Trees.boundingsphereofspheres(SVector(0.0, 0.0), 5.0, SVector(1.0, 0.0), 0.5)
    @test rr ≈ 5.0
    @test norm(cc .- SVector(0.0, 0.0)) <= 1e-12
    # mixed integer/float radii and centers must promote rather than MethodError
    cm, rm = H2Trees.boundingsphereofspheres(SVector(0.0, 0.0), 1, SVector(4, 0.0), 2.0)
    @test norm(cm .- SBT.center(ref)) <= 1e-12
    @test abs(rm - SBT.radius(ref)) <= 1e-12
end

@testset "type stability of adapter return" begin
    rng = MersenneTwister(6)
    points = [SVector(randn(rng), randn(rng), randn(rng)) for _ in 1:16]
    tree = buildtree(
        points;
        builder=BoundingBallTreeBuilder(; splitter=median_split, numsplits=2, minvalues=1),
    )
    node = first(
        Iterators.filter(n -> !H2Trees.isleaf(tree, n), H2Trees.DepthFirstIterator(tree))
    )
    c, r = @inferred H2Trees.boundingsphere(tree, node)
    @test c isa SVector{3,Float64}
    @test r isa Float64
end

@testset "manual 8-child node uses exact SEBB" begin
    # A hand-built one-root/eight-leaf-children tree, independent of any splitter, so the
    # child fanout is exactly 8 (not whatever a binary median split happens to produce) and
    # the true answer is analytically known (plan Section 21 tangent construction).
    sv(x, y, z) = SVector{3,Float64}(x, y, z)

    cstar = sv(1.0, 2.0, -1.0)
    Rstar = 6.0
    dirs = [
        normalize(sv(1, 1, 1)),
        normalize(sv(1, -1, -1)),
        normalize(sv(-1, 1, -1)),
        normalize(sv(-1, -1, 1)),
    ]
    supportradii = [0.2, 0.4, 0.6, 0.8]
    support = [
        SBT.Ball(cstar - (Rstar - supportradii[i]) * dirs[i], supportradii[i]) for
        i in eachindex(dirs)
    ]
    interior = [
        SBT.Ball(cstar + sv(0.1, 0.2, 0.0), 0.5),
        SBT.Ball(cstar + sv(-0.3, 0.1, 0.2), 0.7),
        SBT.Ball(cstar + sv(0.2, -0.4, 0.1), 0.4),
        SBT.Ball(cstar, 1.0),
    ]
    for b in interior
        @test norm(SBT.center(b) - cstar) + SBT.radius(b) < Rstar - 1e-6
    end
    balls = vcat(support, interior)
    @test length(balls) == 8

    # Build the tree: node 1 is the root (firstchild = 2), nodes 2:9 are its leaf children
    # chained by nextsibling, each holding one of `balls` above.
    placeholder = sv(0.0, 0.0, 0.0)
    # Not `H2Trees.Node[...]`: that literal syntax fixes the eltype to the abstract `Node`
    # UnionAll, but `BoundingBallTree`'s constructor requires a concretely-typed
    # `Vector{Node{D}}`. A plain `[...]` literal lets Julia infer the concrete
    # `Node{BoundingBallData{3,Float64}}` element type instead.
    nodes = [H2Trees.Node(H2Trees.BoundingBallData(Int[], placeholder, 0.0, 1), 0, 0, 2)]
    for i in 1:8
        nodeid = i + 1
        nextsib = i == 8 ? 0 : nodeid + 1
        leafdata = H2Trees.BoundingBallData(
            [i], SBT.center(balls[i]), SBT.radius(balls[i]), 2
        )
        push!(nodes, H2Trees.Node(leafdata, nextsib, 1, 0))
    end

    tree = H2Trees.BoundingBallTree(nodes, 1, placeholder, 0.0, [[1], collect(2:9)])
    for i in 1:8
        @test H2Trees.isleaf(tree, i + 1)
    end
    @test !H2Trees.isleaf(tree, 1)

    H2Trees.updateradii!(tree; update=H2Trees.boundingsphere)

    @test norm(H2Trees.center(tree, 1) - cstar) <= 1e-9
    @test abs(H2Trees.radius(tree, 1) - Rstar) <= 1e-9

    # Adapter must match direct SEBB on the same 8 balls, and empty-list/sizehint assumptions
    # about a fixed "typical" fanout must not silently break with exactly 8 children.
    direct = SBT.smallest_enclosing_ball(balls)
    @test norm(H2Trees.center(tree, 1) - SBT.center(direct)) <= 1e-12
    @test abs(H2Trees.radius(tree, 1) - SBT.radius(direct)) <= 1e-12

    # Leaves are untouched by updateradii! (their stored ball is definitionally correct).
    for i in 1:8
        @test H2Trees.center(tree, i + 1) == SBT.center(balls[i])
        @test H2Trees.radius(tree, i + 1) == SBT.radius(balls[i])
    end
end

@testset "multi-way ball tree at large coordinates (regression)" begin
    # A `BoundingBallTree` whose nodes have more than two children is the only shape that
    # exercises `_candidate_many`, and `boundingsphereofspheres`/`median_split` above never do:
    # two children are handled by the exact closed-form two-ball case.
    #
    # For coordinates larger than roughly 1e4 the tangency quadratic has a constant term of
    # order R^2 against a quadratic coefficient of order one, which used to be misread as a
    # degenerate equation. Every candidate for a three- or four-ball support was then discarded
    # and the build died with "SEBB failed to find a numerically valid enclosing ball" -- the
    # failure first seen on the Windows CI runner. It reproduces on any platform once the tree
    # is both wide enough and large enough, which is what this builds.
    for scale in (1.0, 1.0e4, 1.0e6, 1.0e8), seed in 1:3
        rng = MersenneTwister(seed)
        points = [scale * SVector(randn(rng), randn(rng), randn(rng)) for _ in 1:200]
        tree = buildtree(
            points;
            builder=BoundingBallTreeBuilder(;
                splitter=nway_split, numsplits=4, minvalues=8
            ),
        )

        sawmultichild = false
        for node in H2Trees.DepthFirstIterator(tree)
            H2Trees.isleaf(tree, node) && continue
            nchildren = length(collect(H2Trees.children(tree, node)))
            nchildren > 2 && (sawmultichild = true)
            cc = H2Trees.center(tree, node)
            cr = H2Trees.radius(tree, node)
            for child in H2Trees.children(tree, node)
                res = tree_residual(
                    cc, cr, H2Trees.center(tree, child), H2Trees.radius(tree, child)
                )
                @test res <= 1e-9 * max(scale, cr)
            end
        end
        # Without this the loop above could pass vacuously on a tree that never split wide.
        @test sawmultichild
    end
end

end # module TestBoundingBallTreeSEBB
