using Test
using BEAST, CompScienceMeshes
using StaticArrays
using Random
using ParallelKMeans
using H2Trees

@testset "Admissibility predicate" begin
    # The regression this whole predicate exists for.
    #
    # Two independently constructed trees do not share a box grid: a `BlockTree` fits each side's
    # root box to its own space, so the two lattices are offset by a small but genuinely GEOMETRIC
    # amount (measured at ~1.7e-4 on a unit-scale sphere -- about 1e11 times `eps`, so no
    # tolerance tightening could ever have covered it).
    #
    # The legacy centre-distance predicate accepts a pair as near when the centres lie within
    # `sqrt(12)*halfsize`. In 3D that radius is EXACTLY the corner-neighbour centre distance
    # (2h*sqrt(3) == sqrt(12)*h), so corner neighbours sit precisely ON the boundary. Offset them
    # by a hair and they fall outside -- admitted into the FAR field, where the far-field
    # interpolation is then asked to interpolate across two boxes that physically touch, which
    # produces catastrophically wrong matrix entries.
    @testset "offset corner neighbours" begin
        h = 1.0
        a = SVector(0.0, 0.0, 0.0)
        # Diagonal corner neighbour: touches `a` at exactly one corner.
        corner = SVector(2h, 2h, 2h)

        @test H2Trees.boxgap(a, corner, h, h) == 0.0
        @test H2Trees.isneargap(a, corner, h, h, 0)
        # ... and the legacy predicate agrees ONLY in the perfectly aligned case.
        @test H2Trees.isnearhalfsize(a, corner, h, 0)

        # Now offset the second grid the way an independently built tree does.
        for offset in (1e-4, 1e-3, 1e-2)
            shifted = corner .+ offset
            @test H2Trees.boxgap(a, shifted, h, h) > 0          # no longer literally touching
            @test H2Trees.isneargap(a, shifted, h, h, 0)        # but still NEAR -- the fix
            # The legacy predicate loses them: this is the bug, pinned.
            @test !H2Trees.isnearhalfsize(a, shifted, h, 0)
        end
    end

    # The margin is one *smaller* half-size, chosen as the midpoint of the interval that a common
    # grid provably cannot populate: touching neighbours sit at gap 0, and the nearest
    # non-neighbour sits at gap 2h. Pin both endpoints so a future retune cannot silently move the
    # boundary onto a value real geometry takes.
    @testset "margin placement" begin
        h = 1.0
        a = SVector(0.0, 0.0, 0.0)

        # Face-adjacent neighbour: gap exactly 0 -> near.
        @test H2Trees.isneargap(a, SVector(2h, 0.0, 0.0), h, h, 0)
        # Nearest non-neighbour on a shared grid: gap exactly 2h -> far, with a full half-size of
        # clearance above the threshold.
        @test H2Trees.boxgap(a, SVector(4h, 0.0, 0.0), h, h) ≈ 2h
        @test !H2Trees.isneargap(a, SVector(4h, 0.0, 0.0), h, h, 0)
        # Just inside / just outside the margin itself.
        @test H2Trees.isneargap(a, SVector(2h + 0.9h, 0.0, 0.0), h, h, 0)
        @test !H2Trees.isneargap(a, SVector(2h + 1.1h, 0.0, 0.0), h, h, 0)

        # Unequal box sizes: the margin is sized by the SMALLER box, so it means the same thing
        # viewed from either side.
        @test H2Trees.isneargap(a, SVector(3.4, 0.0, 0.0), 1.0, 2.0, 0) ==
            H2Trees.isneargap(SVector(3.4, 0.0, 0.0), a, 2.0, 1.0, 0)

        # `additionalbufferboxes` widens the margin on top of the default one. Probed at gap 1.9,
        # comfortably inside both thresholds rather than sitting on one -- an assertion placed
        # exactly on a boundary is the very failure mode this predicate exists to avoid.
        @test H2Trees.boxgap(a, SVector(3.9h, 0.0, 0.0), h, h) ≈ 1.9h
        @test !H2Trees.isneargap(a, SVector(3.9h, 0.0, 0.0), h, h, 0)
        @test H2Trees.isneargap(a, SVector(3.9h, 0.0, 0.0), h, h, 1)
    end

    # `ballgap`/`nodesize`/`nodegap` are `checkadmissibility`'s BoundingBallTree/KMeansTree
    # counterparts to `boxgap`/`halfsize` -- pinned directly so a future edit to the box path
    # cannot silently leave the ball path broken (that used to fail with a bare `MethodError`
    # from inside `halfsize(::BoundingBallData)`, since only `radius` exists there).
    @testset "ballgap / nodesize / nodegap" begin
        a = SVector(0.0, 0.0, 0.0)

        # Touching balls: gap exactly 0.
        @test H2Trees.ballgap(a, SVector(3.0, 0.0, 0.0), 1.0, 2.0) == 0.0
        # Overlapping balls: still floored at 0, not negative.
        @test H2Trees.ballgap(a, SVector(2.0, 0.0, 0.0), 1.0, 2.0) == 0.0
        # Separated balls: centre distance minus both radii.
        @test H2Trees.ballgap(a, SVector(5.0, 0.0, 0.0), 1.0, 2.0) ≈ 2.0

        # `nodesize`/`nodegap` dispatch on the tree's own shape.
        tree = H2Trees.buildtree(
            [SVector(0.0, 0.0, 0.0), SVector(2.0, 0.0, 0.0)];
            builder=H2Trees.TwoNTreeBuilder(; minhalfsize=0.0, minvalues=0),
        )
        @test H2Trees.nodesize(tree, H2Trees.root(tree)) ==
            H2Trees.halfsize(tree, H2Trees.root(tree))
        @test H2Trees.nodegap(tree, H2Trees.root(tree), tree, H2Trees.root(tree)) == 0.0

        balltree = H2Trees.BoundingBallTree(SVector(0.0, 0.0, 0.0), 1.0)
        @test H2Trees.nodesize(balltree, H2Trees.root(balltree)) ==
            H2Trees.radius(balltree, H2Trees.root(balltree))
        @test H2Trees.nodegap(
            balltree, H2Trees.root(balltree), balltree, H2Trees.root(balltree)
        ) == 0.0
    end
end

@testset "checkadmissibility on BoundingBallTree/KMeansTree" begin
    # `checkadmissibility` must work for ball-shaped trees too, not only `TwoNTree`/`BlockTree` --
    # dispatched via `nodegap`/`nodesize` on the tree's own `treetrait`, mirroring how `isnear`
    # itself already dispatches (`isneargap` for `isTwoNTree`, `isnearradius` for
    # `isBoundingBallTree`).
    Random.seed!(3)
    points = [SVector(randn(), randn(), randn()) for _ in 1:150]

    newsplit(pts, globalpointids, level, numsplits) = (
        [globalpointids[1:2:end], globalpointids[2:2:end]],
        [
            sum(pts[globalpointids[1:2:end]]) / max(1, length(globalpointids[1:2:end])),
            sum(pts[globalpointids[2:2:end]]) / max(1, length(globalpointids[2:2:end])),
        ],
        [1.0, 1.0],
    )
    balltree = H2Trees.buildtree(
        points;
        builder=H2Trees.BoundingBallTreeBuilder(;
            splitter=newsplit, numsplits=2, minvalues=5
        ),
    )
    ballplans = H2Trees.buildplans(balltree)
    ballreport = H2Trees.checkadmissibility(balltree, ballplans; throw=false)
    # `ok` isn't asserted true here: `newsplit` above deliberately ignores spatial locality (it
    # splits by point PARITY, not position), so marginal-gap warnings are expected. The point of
    # this test is that it runs to completion at all -- no `MethodError` from a box-only code path.
    @test ballreport isa H2Trees.AdmissibilityReport

    kmeanstree = H2Trees.KMeansTree(
        points; builder=H2Trees.KMeansTreeBuilder(; numberofclusters=4, minvalues=5)
    )
    kmeansplans = H2Trees.buildplans(kmeanstree)
    kmeansreport = H2Trees.checkadmissibility(kmeanstree, kmeansplans; throw=false)
    @test kmeansreport isa H2Trees.AdmissibilityReport
    @test kmeansreport.ok
    @test isempty(filter(f -> f.severity === :error, kmeansreport.findings))
end

@testset "checkadmissibility" begin
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere6.in")
    )

    @testset "Galerkin plans are admissible" begin
        X = raviartthomas(m)
        tree = H2Trees.buildtree(
            X; builder=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0)
        )
        plans = H2Trees.buildplans(tree)

        report = H2Trees.checkadmissibility(tree, plans; throw=false)
        @test report.ok
        @test isempty(report.findings)
        # `throw=true` must be a no-op on a clean scheme.
        @test H2Trees.checkadmissibility(tree, plans; throw=true) isa
            H2Trees.AdmissibilityReport
    end

    @testset "safe to call nested inside, or concurrently with, another threaded region" begin
        # `checkadmissibility` used to schedule its coverage pass with `Threads.@threads
        # :static`, which errors ("`@threads :static` cannot be used concurrently or nested")
        # if it is ever invoked from inside the caller's own threaded loop, or if two
        # `checkadmissibility` calls happen to run at the same time via separate `@spawn`ed
        # tasks -- both realistic ways to call a public diagnostic function, not just a
        # theoretical concern. Fixed by chunking the coverage pass manually (each chunk's task
        # owns its own scratch buffer) instead of relying on `:static` thread-to-task ownership.
        X = raviartthomas(m)
        tree = H2Trees.buildtree(
            X; builder=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0)
        )
        plans = H2Trees.buildplans(tree)
        expected = H2Trees.checkadmissibility(tree, plans; throw=false)
        @test expected.ok

        nested = fill(false, 4)
        Threads.@threads for i in eachindex(nested)
            nested[i] = H2Trees.checkadmissibility(tree, plans; throw=false).ok
        end
        @test all(nested)

        concurrent = [
            Threads.@spawn(H2Trees.checkadmissibility(tree, plans; throw=false)) for
            _ in 1:4
        ]
        @test all(t -> fetch(t).ok, concurrent)
    end

    @testset "Petrov plans are admissible" begin
        # Genuinely distinct test/trial spaces, so the two BlockTree sides are built independently
        # -- the configuration the corner-neighbour bug lived in.
        X = lagrangec0d1(m)
        Y = lagrangecxd0(m)
        blocktree = H2Trees.buildtree(
            X,
            Y;
            builder=H2Trees.BlockTreeBuilder(;
                test=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0),
                trial=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0),
            ),
        )
        plans = H2Trees.buildplans(blocktree)

        report = H2Trees.checkadmissibility(blocktree, plans; throw=false)
        @test report.ok
        @test isempty(report.findings)
    end

    @testset "detects a plan that contradicts its predicate" begin
        X = lagrangec0d1(m)
        Y = lagrangecxd0(m)
        blocktree = H2Trees.buildtree(
            X,
            Y;
            builder=H2Trees.BlockTreeBuilder(;
                test=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0),
                trial=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0),
            ),
        )
        plans = H2Trees.buildplans(blocktree)

        # A predicate calling everything near turns every scheduled translation into an error --
        # this is the shape of the real bug (pairs in the plan that the predicate says are near),
        # just induced deliberately.
        alwaysnear = tree -> ((a, b, c, d) -> true)
        report = H2Trees.checkadmissibility(
            blocktree, plans; isnear=alwaysnear, coverage=false, throw=false
        )
        @test !report.ok
        @test !isempty(report.findings)
        @test all(f -> f.kind === :nearpairtranslated, report.findings)
        @test all(f -> f.severity === :error, report.findings)
        # Gaps are reported so a caller can see how marginal the offenders were.
        @test all(f -> !ismissing(f.gap), report.findings)

        @test_throws ErrorException H2Trees.checkadmissibility(
            blocktree, plans; isnear=alwaysnear, coverage=false, throw=true
        )
    end

    @testset "marginal far pairs warn without failing" begin
        X = raviartthomas(m)
        tree = H2Trees.buildtree(
            X; builder=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0)
        )
        plans = H2Trees.buildplans(tree)

        # On a shared grid the closest far pairs sit at gap == 2*halfsize exactly, so demanding
        # more than that flags them -- as warnings, which must NOT clear `ok`.
        report = H2Trees.checkadmissibility(
            tree, plans; mingapboxes=3.0, coverage=false, throw=false
        )
        @test report.ok
        @test !isempty(report.findings)
        @test all(f -> f.kind === :marginalgap, report.findings)
        @test all(f -> f.severity === :warning, report.findings)
        # The common-grid minimum: no far pair is closer than two half-sizes.
        @test all(f -> f.gapoverminnodesize >= 2.0 - 1e-9, report.findings)
    end

    @testset "findings are capped" begin
        X = lagrangec0d1(m)
        Y = lagrangecxd0(m)
        blocktree = H2Trees.buildtree(
            X,
            Y;
            builder=H2Trees.BlockTreeBuilder(;
                test=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0),
                trial=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0),
            ),
        )
        plans = H2Trees.buildplans(blocktree)
        alwaysnear = tree -> ((a, b, c, d) -> true)
        report = H2Trees.checkadmissibility(
            blocktree, plans; isnear=alwaysnear, coverage=false, throw=false, maxfindings=7
        )
        @test length(report.findings) <= 7
    end

    # These two negatives exist because the side/level resolution is the subtlest part of the
    # checker: `DisaggregateTranslatePlan` owns its RECEIVING side while `AggregateTranslatePlan`
    # owns its TRANSLATING side, and getting that inversion wrong makes every node in a Petrov
    # plan look mis-levelled. Corrupting a plan deliberately is the only way to prove the check
    # actually fires rather than being vacuously satisfied.
    @testset "detects malformed plans" begin
        X = raviartthomas(m)

        @testset "node that does not exist in its tree" begin
            tree = H2Trees.buildtree(
                X; builder=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0)
            )
            plans = H2Trees.buildplans(tree)
            plan = plans.testdisaggregationplan
            leveldict = first(filter(!isempty, H2Trees.translatingnodes(plan)))
            leveldict[first(keys(leveldict))] = [typemax(Int) - 1]

            report = H2Trees.checkadmissibility(tree, plans; coverage=false, throw=false)
            @test !report.ok
            @test any(f -> f.kind === :unknownnode, report.findings)
        end

        @testset "node filed at the wrong level" begin
            tree = H2Trees.buildtree(
                X; builder=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0)
            )
            plans = H2Trees.buildplans(tree)
            plan = plans.testdisaggregationplan
            # Splice a node from a genuinely different level into a translating list.
            deeplevel = last(H2Trees.levels(plan))
            shallowlevel = first(H2Trees.levels(plan))
            @test deeplevel != shallowlevel
            shallownode = first(H2Trees.disaggregationnodes(plan, shallowlevel))
            deepreceivers = H2Trees.receivingnodes(plan, deeplevel)
            leveldict = H2Trees.translatingnodes(plan)[H2Trees.leveltolevelid(
                plan, deeplevel
            )]
            leveldict[first(deepreceivers)] = [shallownode]

            report = H2Trees.checkadmissibility(tree, plans; coverage=false, throw=false)
            @test !report.ok
            @test any(f -> f.kind === :wronglevel, report.findings)
        end
    end

    # Coverage is read from the PLAN, so it must catch a plan whose far set is wrong even when the
    # predicate is perfectly fine -- an omitted or duplicated translating node. Without these two,
    # the coverage check could be silently vacuous (it would pass just as happily if it were still
    # recomputing both halves from the iterators, which is what it used to do).
    @testset "detects wrong far coverage in the plan" begin
        X = raviartthomas(m)

        function galerkinplansfor(tree)
            return H2Trees.buildplans(tree)
        end

        # Find a (leveldict, receivingnode) that actually schedules translations.
        function firstpopulated(plan)
            for leveldict in H2Trees.translatingnodes(plan)
                for (receivingnode, translating) in leveldict
                    isempty(translating) || return (leveldict, receivingnode)
                end
            end
            return nothing
        end

        @testset "omitted translating node" begin
            tree = H2Trees.buildtree(
                X; builder=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0)
            )
            plans = galerkinplansfor(tree)
            leveldict, receivingnode = firstpopulated(plans.testdisaggregationplan)
            leveldict[receivingnode] = leveldict[receivingnode][2:end]

            report = H2Trees.checkadmissibility(tree, plans; throw=false)
            @test !report.ok
            @test any(f -> f.kind === :coveragegap, report.findings)
        end

        @testset "duplicated translating node" begin
            tree = H2Trees.buildtree(
                X; builder=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0)
            )
            plans = galerkinplansfor(tree)
            leveldict, receivingnode = firstpopulated(plans.testdisaggregationplan)
            leveldict[receivingnode] = [
                leveldict[receivingnode]
                first(leveldict[receivingnode])
            ]

            report = H2Trees.checkadmissibility(tree, plans; throw=false)
            @test !report.ok
            @test any(f -> f.kind === :coverageduplicate, report.findings)
        end
    end

    # `isnear` is diagnostic API, so it takes whichever shape the caller happens to have: a factory
    # `tree -> predicate` (what `H2Trees.isnear()` returns), or the resolved node-comparison
    # callable (what `H2Trees.isnear(A)` returns from a map, and what anyone hand-writing a closure
    # would produce). Both must reach the same verdict.
    @testset "accepts either isnear shape" begin
        @testset "single tree" begin
            X = raviartthomas(m)
            tree = H2Trees.buildtree(
                X; builder=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0)
            )
            plans = H2Trees.buildplans(tree)

            factoryform = H2Trees.checkadmissibility(
                tree, plans; coverage=false, throw=false
            )
            # Resolved, hand-written, three-argument: exactly what a user would type.
            rawclosure =
                (t, nodea, nodeb) -> H2Trees.isnear(t, nodea, nodeb, H2Trees.treetrait(t))
            closureform = H2Trees.checkadmissibility(
                tree, plans; isnear=rawclosure, coverage=false, throw=false
            )

            @test factoryform.ok == closureform.ok
            @test length(factoryform.findings) == length(closureform.findings)
        end

        @testset "block tree" begin
            X = lagrangec0d1(m)
            Y = lagrangecxd0(m)
            blocktree = H2Trees.buildtree(
                X,
                Y;
                builder=H2Trees.BlockTreeBuilder(;
                    test=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0),
                    trial=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0),
                ),
            )
            plans = H2Trees.buildplans(blocktree)

            factoryform = H2Trees.checkadmissibility(
                blocktree, plans; coverage=false, throw=false
            )
            # Resolved, hand-written, four-argument.
            rawclosure =
                (tt, rt, testnode, trialnode) -> H2Trees.isnear(
                    tt,
                    rt,
                    testnode,
                    trialnode,
                    H2Trees.treetrait(tt),
                    H2Trees.treetrait(rt),
                )
            closureform = H2Trees.checkadmissibility(
                blocktree, plans; isnear=rawclosure, coverage=false, throw=false
            )

            @test factoryform.ok == closureform.ok
            @test length(factoryform.findings) == length(closureform.findings)
        end

        @testset "a predicate of neither shape is rejected up front" begin
            X = raviartthomas(m)
            tree = H2Trees.buildtree(
                X; builder=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0)
            )
            plans = H2Trees.buildplans(tree)
            # Two arguments: neither a tree-factory nor a valid node-comparison signature. This must
            # fail immediately with an actionable message, not as a MethodError from inside the
            # translation-pair loop.
            @test_throws ErrorException H2Trees.checkadmissibility(
                tree, plans; isnear=(a, b) -> true, coverage=false, throw=false
            )
        end
    end

    @testset "rejects a tree the plans were not built from" begin
        X = lagrangec0d1(m)
        Y = lagrangecxd0(m)
        blocktree = H2Trees.buildtree(
            X,
            Y;
            builder=H2Trees.BlockTreeBuilder(;
                test=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0),
                trial=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0),
            ),
        )
        plans = H2Trees.buildplans(blocktree)
        # A structurally identical but distinct BlockTree: side resolution is by `===` identity,
        # so this must be refused rather than silently checked against the wrong side.
        other = H2Trees.buildtree(
            X,
            Y;
            builder=H2Trees.BlockTreeBuilder(;
                test=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0),
                trial=H2Trees.TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0),
            ),
        )
        @test_throws ErrorException H2Trees.checkadmissibility(
            other, plans; coverage=false, throw=false
        )
    end
end
