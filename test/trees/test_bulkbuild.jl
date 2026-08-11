using Test
using StaticArrays
using CompScienceMeshes, BEAST
using H2Trees

include(joinpath(pkgdir(H2Trees), "test", "testutils.jl"))

# `bulkbuildtwontree` is the only production path for box trees. These tests pin both internal
# regimes: adaptive `minvalues`/protrusion construction and uniform-depth construction.

function _detpoints(N::Int, n::Int)
    golden = (sqrt(5.0) - 1.0) / 2.0
    return [SVector(ntuple(d -> mod(i * golden^d, 1.0) * 10.0, N)) for i in 1:n]
end

# Test functor that makes the root conformance gate fire.
struct _AlwaysProtrudesFunctor end
(::_AlwaysProtrudesFunctor)(center, halfsize, value) = value == 1 ? 10.0 : 0.0

# Legacy-only protrusion functor; construction must reject this call shape.
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
        # protrusion check; the next testset makes it fire.
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
        # Bounded purely by minhalfsize; every leaf's halfsize must satisfy the bound (or be the
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
        # the general adaptive path, see "coincident points still throw" below); they just share
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
        # `halfsize <= minhalfsize` never becomes true: unbounded recursion instead of a clear
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
        # `bulkbuildtwontree` always calls `protrusion.compute` as `f(center, halfsize, value)`;
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

    @testset "coincident points throw only when they prevent termination" begin
        # Duplicate input points are a caller error and are NOT validated (see
        # `bulkbuildtwontree`'s docstring). The `ArgumentError` fires only where duplicates would
        # actually stop the build terminating: `_uniformseparationdepth` now scans with
        # `minvalues` as its cutoff, so a cell needs MORE than `minvalues` coincident points
        # before subdivision can no longer bring it under the stopping rule.
        duplicate = SVector(0.0, 0.0, 0.0)
        minvalues = 5

        # More duplicates than `minvalues`: no amount of splitting shrinks the cell below the
        # stopping rule, `halfsize` reaches 0 with them still together -> must throw, otherwise
        # construction would recurse forever.
        unresolvable = fill(duplicate, minvalues + 3)
        rootcenter, rootsize = H2Trees.boundingbox(unresolvable)
        rootcenter = SVector(rootcenter...)
        @test_throws ArgumentError H2Trees.bulkbuildtwontree(
            unresolvable,
            rootcenter,
            rootsize,
            zero(rootsize),
            1,
            1,
            minvalues,
            H2Trees.NoProtrusionCheck(),
        )

        # Fewer duplicates than `minvalues`: the cell stops on the `minvalues` rule before the
        # geometry degenerates, so the build succeeds and they share a leaf, matching what the
        # `minvalues == 0` uniform path has always done with coincident points.
        benign = [duplicate, duplicate, SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0)]
        rootcenter, rootsize = H2Trees.boundingbox(benign)
        rootcenter = SVector(rootcenter...)
        tree = H2Trees.bulkbuildtwontree(
            benign,
            rootcenter,
            rootsize,
            zero(rootsize),
            1,
            1,
            minvalues,
            H2Trees.NoProtrusionCheck(),
        )
        @test sort(H2Trees.values(tree, H2Trees.root(tree))) == collect(eachindex(benign))
        # Both copies land in the same leaf, and every input index is accounted for exactly once.
        leafof = Dict(
            value => leaf for leaf in H2Trees.leaves(tree) for
            value in H2Trees.values(H2Trees.data(tree, leaf))
        )
        @test leafof[1] == leafof[2]
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

# Both of these were review findings on the pre-publication `main`.
@testset "root honours the same stop rules as every other node" begin
    # `minvalues` was only consulted inside `_bulkbuildnode!`, i.e. AFTER the root had already
    # been split, so a tree small enough to stop at the root was still subdivided once.
    points = [SVector(float(i), 0.0, 0.0) for i in 1:10]
    tree = H2Trees.buildtree(
        points; builder=H2Trees.TwoNTreeBuilder(; minvalues=70, minhalfsize=0.0)
    )
    @test H2Trees.isleaf(tree, H2Trees.root(tree))
    @test length(collect(H2Trees.DepthFirstIterator(tree))) == 1
    @test sort(H2Trees.values(tree, H2Trees.root(tree))) == collect(eachindex(points))

    # exactly at the boundary the root still stops; one more point and it splits
    atlimit = [SVector(float(i), 0.0, 0.0) for i in 1:8]
    @test H2Trees.isleaf(
        H2Trees.buildtree(
            atlimit; builder=H2Trees.TwoNTreeBuilder(; minvalues=8, minhalfsize=0.0)
        ),
        1,
    )
    @test !H2Trees.isleaf(
        H2Trees.buildtree(
            atlimit; builder=H2Trees.TwoNTreeBuilder(; minvalues=7, minhalfsize=0.0)
        ),
        1,
    )

    # and the root stop must not swallow the protrusion/minhalfsize regimes
    fine = H2Trees.buildtree(
        points; builder=H2Trees.TwoNTreeBuilder(; minvalues=2, minhalfsize=0.0)
    )
    @test !H2Trees.isleaf(fine, H2Trees.root(fine))
end

@testset "unsupported dimensions are rejected, not crashed into" begin
    # `_separationdepth!` counts sectors in a fixed `MVector{8}` under `@inbounds`, which is
    # only wide enough for N <= 3. A 4D input used to index past it and SEGFAULT the process
    # before `HilbertOrdering`'s own dimension check could ever run.
    points4d = [SVector(float(i), float(i % 3), float(i % 5), float(i % 7)) for i in 1:40]
    @test_throws ArgumentError H2Trees.buildtree(
        points4d; builder=H2Trees.TwoNTreeBuilder(; minvalues=10)
    )

    # the supported dimensions keep working
    for N in 1:3
        pts = [SVector{N,Float64}(ntuple(d -> float(i + d), N)) for i in 1:40]
        @test H2Trees.buildtree(pts; builder=H2Trees.TwoNTreeBuilder(; minvalues=10)) isa
            H2Trees.TwoNTree
    end
end
