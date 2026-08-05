using Test
using StaticArrays
using CompScienceMeshes, BEAST
using H2Trees

include(joinpath(pkgdir(H2Trees), "test", "testutils.jl"))

# `bulkbuildtwontree` is the ONLY production construction path for `TwoNTree`/`BlockTree`/
# `SimpleHybridTree` -- it has two internal regimes (see its own docstring): a general adaptive
# path (minvalues/protrusion-driven conformance, bounded by `_uniformseparationdepth` as a safety
# cap, which throws for coincident points) and a uniform-depth bulk path for `minvalues=0`+
# `NoProtrusionCheck` (bounded purely by `minhalfsize`, tolerant of coincident points). These tests
# check properties of each regime directly, rather than comparing against a second
# implementation -- there is no longer a separate point-by-point construction algorithm in
# production to compare against (see git history around the removal of `_addpointbasedtwontree`
# for the equivalence testing that validated this consolidation).

function _detpoints(N::Int, n::Int)
    golden = (sqrt(5.0) - 1.0) / 2.0
    return [SVector(ntuple(d -> mod(i * golden^d, 1.0) * 10.0, N)) for i in 1:n]
end

# `ComputeProtrusionFunctor`'s base implementation always returns `zero(T)`, so sweeping
# `maxprotrusion` against it never actually trips the protrusion check -- see the
# "protrusion actually trips" testset below, which uses this functor instead to confirm the
# root-blocks-everyone gate genuinely fires (not just executes).
struct _AlwaysProtrudesFunctor end
(::_AlwaysProtrudesFunctor)(center, halfsize, value) = value == 1 ? 10.0 : 0.0

# Deliberately implements ONLY the legacy `f(tree, node, value)` call form -- no longer a
# supported construction-time protrusion interface (see `bulkbuildtwontree`'s docstring), so using
# it must fail loudly rather than being silently accepted or falling back to something else.
struct _LegacyOnlyProtrusionFunctor end
function (::_LegacyOnlyProtrusionFunctor)(tree, node::Int, value::Int)
    return 0.1 * value / max(1, H2Trees.numberofvalues(tree, node))
end

function _valuecoverage(tree)
    return sort!(reduce(vcat, H2Trees.values.(Ref(tree), H2Trees.leaves(tree))))
end

@testset "bulkbuildtwontree" begin
    @testset "general adaptive path: synthetic point clouds, minvalues-only" begin
        for N in (1, 2, 3)
            points = _detpoints(N, 90)
            rootcenter, rootsize = H2Trees.boundingbox(points)
            rootcenter = SVector(rootcenter...)
            for minvalues in (1, 3, 10, 20)
                for minlevel in (1, 2)
                    for root in (1, 2)
                        tree = H2Trees.bulkbuildtwontree(
                            points,
                            rootcenter,
                            rootsize,
                            zero(rootsize),
                            minlevel,
                            root,
                            minvalues,
                            H2Trees.NoProtrusionCheck(),
                        )
                        @test _valuecoverage(tree) == 1:length(points)
                        @test H2Trees.root(tree) == root
                        @test H2Trees.levels(tree)[1] == minlevel
                        for leaf in H2Trees.leaves(tree)
                            leaf == H2Trees.root(tree) && continue
                            @test H2Trees.numberofvalues(tree, H2Trees.parent(tree, leaf)) >
                                minvalues
                        end
                    end
                end
            end
        end
    end

    @testset "TreeIndex/near-far/plan validity on a bulk-built tree" begin
        # Beyond value coverage, the bulk-built tree must be internally valid on its own terms:
        # TreeIndex bookkeeping, near/far well-separatedness, and plan/checkadmissibility all have
        # to work against it exactly as they do for any other tree.
        points = _detpoints(3, 120)
        rootcenter, rootsize = H2Trees.boundingbox(points)
        rootcenter = SVector(rootcenter...)
        tree = H2Trees.buildtree(
            points;
            builder=H2Trees.TwoNTreeBuilder(; minhalfsize=rootsize / 16, minvalues=1),
        )

        @test H2Trees.nodesatlevel(tree) == H2Trees.treeindex(tree).nodes_by_level
        @test H2Trees.treeindex(tree).leaves == H2Trees.leaves(tree)
        @test H2Trees.depthfirstnodes(tree) ==
            collect(Int, H2Trees.DepthFirstIterator(tree))
        @test _valuecoverage(tree) == 1:length(points)
        @test TestingUtils.testwellseparatedness(tree)

        plans = H2Trees.buildplans(tree; builder=H2Trees.PlanBuilder())
        report = H2Trees.checkadmissibility(tree, plans; throw=false)
        @test report.ok
    end

    @testset "general adaptive path: minhalfsize interacting with minvalues" begin
        points = _detpoints(3, 200)
        rootcenter, rawrootsize = H2Trees.boundingbox(points)
        rootcenter = SVector(rootcenter...)
        for minhalfsize in (0.0, 0.05, 0.2)
            for minvalues in (1, 5, 20)
                adjustedhalfsize = H2Trees.roothalfsize(rawrootsize, minhalfsize)
                tree = H2Trees.bulkbuildtwontree(
                    points,
                    rootcenter,
                    adjustedhalfsize,
                    minhalfsize,
                    1,
                    1,
                    minvalues,
                    H2Trees.NoProtrusionCheck(),
                )
                @test _valuecoverage(tree) == 1:length(points)
                for leaf in H2Trees.leaves(tree)
                    @test H2Trees.halfsize(tree, leaf) >= minhalfsize
                end
            end
        end
    end

    @testset "protrusion enabled but structurally inert (ComputeProtrusionFunctor's base case)" begin
        # `ComputeProtrusionFunctor()` always returns 0, so this sweep never actually trips the
        # protrusion check -- see the next testset for the case where it genuinely fires.
        points = _detpoints(3, 90)
        rootcenter, rootsize = H2Trees.boundingbox(points)
        rootcenter = SVector(rootcenter...)
        for maxprotrusion in (0.01, 0.3, 5.0)
            for minvalues in (0, 5, 15)
                protrusion = H2Trees.ProtrusionCheck(;
                    max=maxprotrusion, compute=H2Trees.ComputeProtrusionFunctor()
                )
                tree = H2Trees.bulkbuildtwontree(
                    points,
                    rootcenter,
                    rootsize,
                    zero(rootsize),
                    1,
                    1,
                    minvalues,
                    protrusion,
                )
                @test _valuecoverage(tree) == 1:length(points)
            end
        end
    end

    @testset "protrusion actually trips: the root-blocks-everyone gate fires for real" begin
        points = _detpoints(3, 90)
        rootcenter, rootsize = H2Trees.boundingbox(points)
        rootcenter = SVector(rootcenter...)
        protrusion = H2Trees.ProtrusionCheck(; max=1.0, compute=_AlwaysProtrudesFunctor())

        tree = H2Trees.bulkbuildtwontree(
            points, rootcenter, rootsize, zero(rootsize), 1, 1, 5, protrusion
        )
        @test H2Trees.numberofnodes(tree) == 1
        @test _valuecoverage(tree) == 1:length(points)
    end

    @testset "uniform-depth bulk path: minvalues=0 + NoProtrusionCheck" begin
        # Bounded purely by minhalfsize -- every leaf's halfsize must satisfy the bound (or be the
        # root, if the root itself is already small enough), and value coverage must hold, without
        # ever computing `_uniformseparationdepth`.
        points = _detpoints(3, 150)
        rootcenter, rootsize = H2Trees.boundingbox(points)
        rootcenter = SVector(rootcenter...)
        for minhalfsize in (rootsize / 4, rootsize / 16, rootsize / 64)
            tree = H2Trees.bulkbuildtwontree(
                points,
                rootcenter,
                rootsize,
                minhalfsize,
                1,
                1,
                0,
                H2Trees.NoProtrusionCheck(),
            )
            @test _valuecoverage(tree) == 1:length(points)
            for leaf in H2Trees.leaves(tree)
                @test H2Trees.halfsize(tree, leaf) <= minhalfsize ||
                    leaf == H2Trees.root(tree)
            end
        end

        # The whole point of this regime: coincident points at minhalfsize=0 don't throw (unlike
        # the general adaptive path, see "coincident points still throw" below) -- they just share
        # a leaf once halfsize underflows to exactly 0.0.
        duplicate = SVector(0.0, 0.0, 0.0)
        coincidentpoints = [duplicate, duplicate, SVector(1.5, -2.0, 0.25)]
        ccenter, csize = H2Trees.boundingbox(coincidentpoints)
        tree = H2Trees.bulkbuildtwontree(
            coincidentpoints,
            SVector(ccenter...),
            csize,
            zero(csize),
            1,
            1,
            0,
            H2Trees.NoProtrusionCheck(),
        )
        @test _valuecoverage(tree) == 1:3
    end

    @testset "invalid protrusion argument is rejected, not silently ignored" begin
        # `protrusion.compute` is public/user-extensible, so `bulkbuildtwontree` must validate its
        # type up front and throw a clear `ArgumentError`, rather than silently treating an
        # invalid `protrusion` as "no protrusion".
        points = _detpoints(3, 50)
        rootcenter, rootsize = H2Trees.boundingbox(points)
        rootcenter = SVector(rootcenter...)
        @test_throws ArgumentError H2Trees.bulkbuildtwontree(
            points, rootcenter, rootsize, zero(rootsize), 1, 1, 5, :notarealprotrusion
        )
        @test_throws ArgumentError H2Trees.buildtree(
            points;
            builder=H2Trees.TwoNTreeBuilder(; minvalues=5, protrusion=:notarealprotrusion),
        )
    end

    @testset "invalid minhalfsize is rejected, not silently accepted" begin
        # In the minvalues=0 + NoProtrusionCheck uniform-bulk regime, minhalfsize is the ONLY
        # thing bounding recursion (`_uniformseparationdepth`'s depth cap is deliberately never
        # computed for this regime). A negative or non-finite minhalfsize would mean
        # `halfsize <= minhalfsize` never becomes true -- unbounded recursion instead of a clear
        # error. The public `TwoNTree(positions; builder)` entry point already rejects a negative
        # minhalfsize indirectly (a confusing `DomainError` from `roothalfsize`'s `log2` call,
        # before `bulkbuildtwontree` is ever reached), but `bulkbuildtwontree`/the direct
        # `TwoNTree(center, points, halfsize, minhalfsize; ...)` constructor bypass that and must
        # validate on their own.
        points = _detpoints(3, 50)
        rootcenter, rootsize = H2Trees.boundingbox(points)
        rootcenter = SVector(rootcenter...)
        for bad in (-0.1, NaN, -Inf, Inf)
            @test_throws ArgumentError H2Trees.bulkbuildtwontree(
                points, rootcenter, rootsize, bad, 1, 1, 0, H2Trees.NoProtrusionCheck()
            )
        end
        # 0.0 (the default) and any finite non-negative value remain valid.
        @test H2Trees.bulkbuildtwontree(
            points, rootcenter, rootsize, 0.0, 1, 1, 0, H2Trees.NoProtrusionCheck()
        ) isa H2Trees.TwoNTree
    end

    @testset "legacy f(tree, node, value) protrusion functors are no longer supported" begin
        # `bulkbuildtwontree` always calls `protrusion.compute` as `f(center, halfsize, value)` --
        # there is no fallback for a functor that only implements the older `f(tree, node, value)`
        # shape; using one must fail (a `MethodError` from the mismatched call), not silently
        # construct a tree via some other mechanism.
        points = _detpoints(3, 90)
        rootcenter, rootsize = H2Trees.boundingbox(points)
        rootcenter = SVector(rootcenter...)
        protrusion = H2Trees.ProtrusionCheck(;
            max=0.3, compute=_LegacyOnlyProtrusionFunctor()
        )
        @test_throws MethodError H2Trees.bulkbuildtwontree(
            points, rootcenter, rootsize, zero(rootsize), 1, 1, 5, protrusion
        )
        @test_throws MethodError H2Trees.buildtree(
            points; builder=H2Trees.TwoNTreeBuilder(; minvalues=5, protrusion=protrusion)
        )
    end

    @testset "coincident points still throw" begin
        duplicate = SVector(0.0, 0.0, 0.0)
        points = [duplicate, duplicate]
        rootcenter, rootsize = H2Trees.boundingbox(points)
        rootcenter = SVector(rootcenter...)
        @test_throws ArgumentError H2Trees.bulkbuildtwontree(
            points,
            rootcenter,
            rootsize,
            zero(rootsize),
            1,
            1,
            5,
            H2Trees.NoProtrusionCheck(),
        )
    end

    @testset "BEAST-space protrusion path" begin
        m = CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere.in")
        )
        X = raviartthomas(m)
        points = BEAST.positions(X)
        rootcenter, rootsize = H2Trees.boundingbox(points)
        rootcenter = SVector(rootcenter...)
        beastprotrusion = H2Trees.BEASTProtrusionFunctor(X)
        for maxprotrusion in (0.1, 0.5, 1.0)
            for minvalues in (0, 10, 100)
                protrusion = H2Trees.ProtrusionCheck(;
                    max=maxprotrusion, compute=beastprotrusion
                )
                tree = H2Trees.bulkbuildtwontree(
                    points,
                    rootcenter,
                    rootsize,
                    zero(rootsize),
                    1,
                    1,
                    minvalues,
                    protrusion,
                )
                @test _valuecoverage(tree) == 1:length(points)
                _maxprotrusion = H2Trees.maxprotrusion(
                    tree; computeprotrusion=beastprotrusion
                )
                @test maximum(_maxprotrusion) <= maxprotrusion
            end
        end
    end
end
