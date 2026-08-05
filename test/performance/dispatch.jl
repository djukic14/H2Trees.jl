# JET runtime-dispatch checks, scoped to H2Trees' own code (`target_modules=(H2Trees,)`) so a
# workload's dependencies (Graphs, Metis, BoundingSphere, ...) don't add unrelated noise.
#
# JET's abstract interpretation reports a small, fixed set of pre-existing false positives that
# show up for essentially any workload going through the affected code -- documented and allowlisted
# below rather than chased to zero, per the plan's own instruction to pin strict checks to one
# Julia version and not treat JET noise as a blocking requirement. Any report NOT matching the
# allowlist fails the test: that is a real, new dispatch issue.

using Test
using JET
using H2Trees
using StaticArrays
using Graphs

# `DepthFirstIterator`'s `iterate` only calls `node(info)` after checking `!isnothing(info)`, but
# that check is a value-level dispatch on a custom `Base.isnothing(::NodeInformation)` override --
# JET's abstract interpreter cannot use it to narrow `info`'s type away from
# `NodeInformation{Nothing}`, so it reports the (never actually reached at runtime) case where
# `node(::NodeInformation{Nothing})` would index into `nothing`. Shows up for every workload that
# touches `DepthFirstIterator` (directly or via `leaves`/`values`/tree construction), which is
# effectively all of them.
_isknownnodeinformationreport(r) = occursin("NodeInformation{Nothing}", sprint(show, r))

# `AggregatePlan`'s `aggregationlevels` field is a `Union` of range types; JET's union-split
# analysis reports the non-matching half of the split as a "no method" error even though the
# actual runtime value always resolves to one concrete member. Pre-existing, not specific to any
# workload added here.
_isknownaggregateplanreport(r) = occursin("AggregatePlan(", sprint(show, r))

# `H2ClusterTree`'s callable syntax (`tree(node::Int)`) makes JET's union-split analysis consider
# (and reject) a branch where the local named `tree` is an `Int64` being called as a functor --
# a shape no real call in `checkadmissibility` ever takes.
_isknowntreecallreport(r) = occursin("(::Int64)(::Int64)", sprint(show, r))

function _isknownreport(r)
    return _isknownnodeinformationreport(r) ||
           _isknownaggregateplanreport(r) ||
           _isknowntreecallreport(r)
end

function perf_unexpectedreports(f, argtypes)
    result = JET.report_call(f, argtypes; target_modules=(H2Trees,))
    return filter(!_isknownreport, JET.get_reports(result))
end

@testset "Construction dispatch" begin
    pts = perf_points(3)
    trial = perf_trial_points(3)
    g = perf_graph(length(pts))
    w = perf_weights(length(pts))
    gf = perf_forest_graph(100, 100)
    wf = perf_weights(200)

    @test isempty(perf_unexpectedreports(perf_buildtwontree, (typeof(pts),)))
    @test isempty(perf_unexpectedreports(perf_buildblocktree, (typeof(pts), typeof(trial))))
    @test isempty(perf_unexpectedreports(perf_buildboundingballtree, (typeof(pts),)))
    @test isempty(perf_unexpectedreports(perf_buildkmeanstree, (typeof(pts),)))
    @test isempty(
        perf_unexpectedreports(perf_buildmetistree, (typeof(pts), typeof(g), typeof(w)))
    )
    @test isempty(
        perf_unexpectedreports(perf_buildmetisforest, (typeof(pts), typeof(gf), typeof(wf)))
    )
    @test isempty(perf_unexpectedreports(perf_buildsimplehybridtree, (typeof(pts),)))
end

@testset "Iterator dispatch" begin
    pts = perf_points(3)
    tree = perf_buildtwontree(pts)
    T = typeof(tree)

    @test isempty(perf_unexpectedreports(perf_depthfirstsum, (T,)))
    @test isempty(perf_unexpectedreports(perf_childrensum, (T,)))
    @test isempty(perf_unexpectedreports(perf_parentupwardssum, (T,)))
    @test isempty(perf_unexpectedreports(perf_levelsum, (T,)))
    @test isempty(perf_unexpectedreports(perf_nearfarsum, (T,)))
end

@testset "Plan / checkadmissibility dispatch" begin
    pts = perf_points(3)
    trial = perf_trial_points(3)
    plantree = perf_buildplantree(pts)
    planblock = perf_buildplanblocktree(pts, trial)
    galerkinplans = perf_buildgalerkinplans(plantree)
    petrovplans = perf_buildpetrovplans(planblock)

    @test isempty(perf_unexpectedreports(perf_buildgalerkinplans, (typeof(plantree),)))
    @test isempty(perf_unexpectedreports(perf_buildpetrovplans, (typeof(planblock),)))
    @test isempty(
        perf_unexpectedreports(
            perf_checkadmissibility, (typeof(plantree), typeof(galerkinplans))
        ),
    )
    @test isempty(
        perf_unexpectedreports(
            perf_checkadmissibility, (typeof(planblock), typeof(petrovplans))
        ),
    )
end
