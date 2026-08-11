module TestNearLists

using Test
using StaticArrays
using H2Trees

# A predicate H2Trees has no `nearcandidatepredicate` method for, so it must fall back.
struct OpaqueNear{N}
    isnear::N
end
(f::OpaqueNear)(tree, nodea::Int, nodeb::Int) = f.isnear(tree, nodea, nodeb)

# A predicate that only NARROWS its inner one, the shape H2Factory's
# GalerkinSymmetricIsNearFunctor has. Narrowing is safe to filter a cache with.
struct NarrowingNear{N}
    isnear::N
end
function (f::NarrowingNear)(tree, nodea::Int, nodeb::Int)
    return f.isnear(tree, nodea, nodeb) && (!H2Trees.isleaf(tree, nodeb) || nodea > nodeb)
end
function H2Trees.nearcandidatepredicate(f::NarrowingNear)
    return H2Trees.nearcandidatepredicate(f.isnear)
end

# Deliberately spatial: median along the widest axis, with real enclosing balls.
function mediansplit(points, globalpointids, level, numsplits)
    pts = points[globalpointids]
    axis = argmax([
        maximum(p -> p[d], pts) - minimum(p -> p[d], pts) for d in 1:length(first(pts))
    ])
    order = sortperm(globalpointids; by=i -> points[i][axis])
    mid = max(1, length(order) ÷ 2)
    groups = filter(
        !isempty, [globalpointids[order[1:mid]], globalpointids[order[(mid + 1):end]]]
    )
    centers = [sum(points[g]) / length(g) for g in groups]
    radii = [
        maximum(i -> H2Trees.norm(points[i] - c), g) for (g, c) in zip(groups, centers)
    ]
    return groups, centers, radii
end

"""
Order-insensitive view of a `nearinteractions` result: the sorted set of `(values, nearvalues)` pairs, which is all that is well defined once the leaf loop is threaded.
"""
pairset((values, nearvalues)) = sort(collect(zip(values, nearvalues)))

function boxpoints(n, seed)
    # deterministic, no Random dependency: a shuffled lattice-ish cloud
    return [
        SVector(sin(seed + 1.7i), cos(seed + 2.3i), sin(seed + 0.9i) * cos(seed + 1.1i)) for
        i in 1:n
    ]
end

@testset "supportsnearlists gates on tree shape" begin
    points = boxpoints(600, 3)
    boxtree = buildtree(points; builder=TwoNTreeBuilder(; minvalues=20))
    balltree = buildtree(
        points;
        builder=BoundingBallTreeBuilder(; splitter=mediansplit, numsplits=2, minvalues=20),
    )
    blocktree = H2Trees.BlockTree(boxtree, boxtree)

    @test H2Trees.supportsnearlists(boxtree)
    @test !H2Trees.supportsnearlists(balltree)
    @test !H2Trees.supportsnearlists(blocktree)

    @test H2Trees.nearlistcache(boxtree, H2Trees.isnear()(boxtree)) isa
        H2Trees.NearListCache
    @test H2Trees.nearlistcache(balltree, H2Trees.isnear()(balltree)) === nothing
    @test H2Trees.nearlistcache(blocktree, H2Trees.isnear()(blocktree)) === nothing

    # An unrecognised predicate falls back even on a supported shape.
    @test H2Trees.nearlistcache(boxtree, OpaqueNear(H2Trees.isnear()(boxtree))) === nothing

    # ... and the plan builder applies the same gate, which is what keeps a ball tree's
    # translations from being silently dropped (see `supportsnearlists`).
    for tree in (balltree, blocktree)
        iterator = H2Trees._TranslatingFunctor(
            tree, H2Trees.WellSeparatedIterator(; isnear=H2Trees.isnear())(tree)
        )
        @test H2Trees._translatinglists(tree, iterator) === nothing
    end
end

# `isnearradius` used to compare against `2max(r1, r2)`, whose effective margin `|r1 - r2|` is
# not monotone up the tree: this pair of balls was near while its parents were far, even with
# the children strictly inside their parents. The SEBB work replaced that with the
# sum-of-radii form, which is monotone under containment. Kept as a regression pin, because it
# is the exact configuration that made the plan sweep drop translations on ball trees.
@testset "ball near predicate is level-monotone under containment" begin
    childa, radiusa = SVector(0.0, 0.0, 0.0), 0.1
    childb, radiusb = SVector(2.0, 0.0, 0.0), 1.0
    parenta, radiusparenta = SVector(-0.9, 0.0, 0.0), 1.0
    parentb, radiusparentb = SVector(2.0, 0.0, 0.0), 1.0

    # each child ball really is inside its parent ball
    @test H2Trees.norm(childa - parenta) + radiusa <= radiusparenta + 1e-12
    @test H2Trees.norm(childb - parentb) + radiusb <= radiusparentb + 1e-12

    # the old `2max` form called the children near (2.0 <= 2*max(0.1,1.0)) and the parents far;
    # the sum-of-radii form calls both far, so the implication holds vacuously here
    @test !H2Trees.isnearradius(childa, childb, radiusa, radiusb)
    @test !H2Trees.isnearradius(parenta, parentb, radiusparenta, radiusparentb)

    # and it holds non-vacuously: bring the children close enough to be near, and the parents
    # must be near too. Under `2max` this is exactly what could fail.
    for scale in (0.2, 0.5, 0.9, 1.0, 1.05)
        nearchildb = SVector(scale * (radiusa + radiusb), 0.0, 0.0)
        # keep the child inside its parent while moving it
        movedparentb = nearchildb
        H2Trees.isnearradius(childa, nearchildb, radiusa, radiusb) || continue
        @test H2Trees.isnearradius(parenta, movedparentb, radiusparenta, radiusparentb)
    end
end

@testset "cached near lists equal the full-level scan" begin
    for (label, points, minvalues) in (
        ("uniform", boxpoints(900, 1), 25),
        # non-uniform, so leaves land on different levels and the ancestor branch runs
        ("clustered", vcat(boxpoints(700, 2) .* 0.02, boxpoints(300, 5)), 15),
    )
        tree = buildtree(points; builder=TwoNTreeBuilder(; minvalues=minvalues))
        for (pname, predicate) in (
            "default" => H2Trees.isnear()(tree),
            "additionalbufferboxes=1" => H2Trees.isnear(; additionalbufferboxes=1)(tree),
            "narrowing wrapper" => NarrowingNear(H2Trees.isnear()(tree)),
        )
            cache = H2Trees.nearlistcache(tree, predicate)
            @test cache isa H2Trees.NearListCache

            # `_LeafPredicateFunctor` narrows the predicate and is NOT level-monotone, so it
            # must filter these lists rather than seed a sweep.
            leafpredicate = H2Trees._LeafPredicateFunctor(predicate)

            for node in H2Trees.DepthFirstIterator(tree)
                for filterpredicate in (predicate, leafpredicate)
                    scanned = collect(
                        Int, H2Trees.NearNodeIterator(tree, node; isnear=filterpredicate)
                    )
                    cached = collect(
                        Int,
                        H2Trees.NearNodeIterator(
                            tree, node; isnear=filterpredicate, nearlists=cache
                        ),
                    )
                    @test sort(scanned) == sort(cached)
                end
            end

            for leaf in H2Trees.leaves(tree)
                @test H2Trees.nearnodevalues(tree, leaf; isnear=predicate) ==
                    H2Trees.nearnodevalues(tree, leaf; isnear=predicate, nearlists=cache)
            end

            # `nearinteractions` fills its outputs from a `@threads` loop under a lock, so
            # which leaf lands at which index is nondeterministic with more than one thread.
            # Only the SET of (values, nearvalues) pairs is defined; comparing the returned
            # vectors directly passes serially and fails at -t4.
            @test pairset(H2Trees.nearinteractions(tree; isnear=predicate)) ==
                pairset(H2Trees.nearinteractions(tree; isnear=OpaqueNear(predicate)))
        end
        @test label isa String
    end
end

end
