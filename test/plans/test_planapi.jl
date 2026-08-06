using Test
using BEAST, CompScienceMeshes
using StaticArrays
using LinearAlgebra
using H2Trees

@testset "Plan level validation" begin
    @test H2Trees._validatedaggregationlevels([4, 3, 2]) == 4:-1:2
    @test H2Trees._validateddisaggregationlevels([2, 3, 4]) == 2:4
    @test_throws(
        ArgumentError(
            "aggregation plan levels must be contiguous and descending, got [4, 2]"
        ),
        H2Trees._validatedaggregationlevels([4, 2])
    )
    @test_throws(
        ArgumentError(
            "disaggregation plan levels must be contiguous and ascending, got [2, 4]"
        ),
        H2Trees._validateddisaggregationlevels([2, 4])
    )
end

@testset "Translate plan construction: negative cases" begin
    # `AggregateTranslatePlan`/`DisaggregateTranslatePlan` share one internal builder
    # (`_buildtranslateplan`); these pin the invariants that builder must preserve.
    testpts = [SVector(0.0 + 0.1i, 0.0, 0.0) for i in 0:29]
    trialpts = [SVector(5.0 + 0.1i, 0.0, 0.0) for i in 0:29]
    block = H2Trees.buildtree(
        testpts,
        trialpts;
        builder=BlockTreeBuilder(;
            test=TwoNTreeBuilder(; minhalfsize=0.0, minvalues=3),
            trial=TwoNTreeBuilder(; minhalfsize=0.0, minvalues=3),
        ),
    )
    tfiterator = H2Trees.TranslatingNodesIterator(; isnear=H2Trees.isnear())

    # A BlockTree must be rejected by the single-tree constructor -- the caller must specify
    # which side via the two-tree form instead.
    @test_throws ArgumentError H2Trees.AggregateTranslatePlan(block, tfiterator(block))
    @test_throws ArgumentError H2Trees.DisaggregateTranslatePlan(block, tfiterator(block))

    # A tree small/coarse enough that everything is near (no possible translating pair) must
    # still refuse to build an `AggregateTranslatePlan` rather than silently returning an empty
    # one.
    tinytree = buildtree(
        [SVector(0.0, 0.0, 0.0), SVector(0.5, 0.0, 0.0), SVector(0.0, 0.5, 0.0)];
        builder=TwoNTreeBuilder(; minhalfsize=0.0, minvalues=10),
    )
    @test_throws(
        ArgumentError("Empty AggregatePlan not supported."),
        H2Trees.AggregateTranslatePlan(tinytree, tfiterator(tinytree))
    )
end

# This testset (and its Petrov counterpart below) deliberately calls the internal
# `_buildgalerkinplans`/`_buildpetrovplans` rather than the public `buildplans`: the whole point
# here is comparing THEIR output field-by-field against plans assembled by hand from the
# low-level constructors (`AggregatePlan`, `AggregateTranslatePlan`, ...), plus custom
# `aggregatenode` predicates (`aggregateallnodes()`/`aggregaterootonly()`) that exercise
# constructor internals `buildplans`/`PlanBuilder` wouldn't otherwise reach in one call. This is
# Bucket C in the review-followup plan's classification -- an internal-contract test, kept
# private on purpose.
@testset "Galerkin Plan" begin
    λ = 1.0
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere6.in")
    )
    X = raviartthomas(m)
    tree = buildtree(X; builder=TwoNTreeBuilder(; minhalfsize=λ / 10, minvalues=0))

    TFIterator = H2Trees.TranslatingNodesIterator(; isnear=H2Trees.isnear())
    aggregatenode = H2Trees.istranslatingnode(; TranslatingNodesIterator=TFIterator)

    trialaggregateplan = H2Trees.AggregatePlan(tree, aggregatenode(tree))
    testdisaggregatetranslateplan = H2Trees.DisaggregateTranslatePlan(
        tree, TFIterator(tree)
    )
    testaggregatetranslateplan = H2Trees.AggregateTranslatePlan(tree, TFIterator(tree))
    trialdisaggregateplan = H2Trees.DisaggregatePlan(tree, aggregatenode(tree))

    aggregateallplans = H2Trees._buildgalerkinplans(
        tree, H2Trees.aggregateallnodes(), TFIterator, H2Trees.AggregateMode()
    )
    @test sum(aggregateallplans.trialaggregationplan.storenode) == length(tree.nodes)

    aggregaterootplans = H2Trees._buildgalerkinplans(
        tree, H2Trees.aggregaterootonly(), TFIterator, H2Trees.AggregateMode()
    )
    @test H2Trees.aggregaterootonly() isa H2Trees.AggregateRootOnlyFunctor
    @test H2Trees.AggregateRootFunctor === H2Trees.AggregateRootOnlyFunctor
    @test H2Trees.aggregaterootonly()(tree)(H2Trees.root(tree))
    @test sum(aggregaterootplans.trialaggregationplan.storenode) == 1
    @test aggregaterootplans.trialaggregationplan.storenode[H2Trees.root(tree)] == 1
    plans = H2Trees._buildgalerkinplans(
        tree, aggregatenode, TFIterator, H2Trees.AggregateMode()
    )

    @test plans isa H2Trees.PlanSet
    @test H2Trees.isgalerkin(plans)
    @test !H2Trees.ispetrov(plans)
    @test H2Trees.tree(plans) === tree
    @test H2Trees.trialaggregationplan(plans) === plans.trialaggregationplan
    @test H2Trees.testdisaggregationplan(plans) === plans.testdisaggregationplan
    @test H2Trees.testaggregationplan(plans) === plans.testaggregationplan
    @test H2Trees.trialdisaggregationplan(plans) === plans.trialdisaggregationplan
    @test H2Trees.relevantlevels(plans) == plans.relevantlevels
    @test hasproperty(plans, :tree)
    @test hasproperty(plans, :family)
    @test H2Trees.planfamily(plans) isa H2Trees.GalerkinPlanFamily
    @test keys(plans) == keys((
        trialaggregationplan=plans.trialaggregationplan,
        testdisaggregationplan=plans.testdisaggregationplan,
        testaggregationplan=plans.testaggregationplan,
        trialdisaggregationplan=plans.trialdisaggregationplan,
        relevantlevels=plans.relevantlevels,
    ))
    @test plans[:trialaggregationplan] === plans.trialaggregationplan
    @test collect(values(plans))[begin] === plans.trialaggregationplan
    @test NamedTuple(plans) == (
        trialaggregationplan=plans.trialaggregationplan,
        testdisaggregationplan=plans.testdisaggregationplan,
        testaggregationplan=plans.testaggregationplan,
        trialdisaggregationplan=plans.trialdisaggregationplan,
        relevantlevels=plans.relevantlevels,
    )
    @test (; plans...) == NamedTuple(plans)
    @test merge((before=:before,), plans) == merge((before=:before,), NamedTuple(plans))
    @test merge(plans, (after=:after,)) == merge(NamedTuple(plans), (after=:after,))
    @test H2Trees.ownedtree(plans.testdisaggregationplan) === tree
    @test H2Trees.receivingtree(plans, plans.testdisaggregationplan) === tree
    @test H2Trees.translatingtree(plans, plans.testdisaggregationplan) === tree
    @test H2Trees.receivingtree(tree, plans.testaggregationplan) === tree
    @test H2Trees.translatingtree(tree, plans.testaggregationplan) === tree

    builtplans = H2Trees.buildplans(
        tree;
        builder=H2Trees.PlanBuilder(;
            aggregationmode=H2Trees.AggregateMode(), isnear=H2Trees.isnear()
        ),
    )
    @test builtplans isa H2Trees.PlanSet
    @test H2Trees.relevantlevels(builtplans) == H2Trees.relevantlevels(plans)

    ptrialaggregateplan = plans.trialaggregationplan
    ptestdisaggregatetranslateplan = plans.testdisaggregationplan
    ptestaggregatetranslateplan = plans.testaggregationplan
    ptrialdisaggregateplan = plans.trialdisaggregationplan

    @test ptrialaggregateplan.levels == trialaggregateplan.levels
    @test ptrialaggregateplan.storenode == trialaggregateplan.storenode
    @test ptrialaggregateplan.tree == trialaggregateplan.tree
    @test ptrialaggregateplan.rootoffset == trialaggregateplan.rootoffset
    for level in trialaggregateplan.levels
        @test sort(H2Trees.aggregationnodes(ptrialaggregateplan, level)) ==
            sort(H2Trees.aggregationnodes(trialaggregateplan, level))
    end

    @test ptestdisaggregatetranslateplan.levels == testdisaggregatetranslateplan.levels
    @test ptestdisaggregatetranslateplan.isdisaggregationnode ==
        testdisaggregatetranslateplan.isdisaggregationnode
    @test ptestdisaggregatetranslateplan.tree == testdisaggregatetranslateplan.tree
    @test ptestdisaggregatetranslateplan.rootoffset ==
        testdisaggregatetranslateplan.rootoffset
    for level in testdisaggregatetranslateplan.levels
        @test sort(H2Trees.disaggregationnodes(ptestdisaggregatetranslateplan, level)) ==
            sort(H2Trees.disaggregationnodes(testdisaggregatetranslateplan, level))
    end
    @test sort(collect(keys(testdisaggregatetranslateplan.translatingnodes))) ==
        sort(collect(keys(ptestdisaggregatetranslateplan.translatingnodes)))
    for i in eachindex(testdisaggregatetranslateplan.translatingnodes)
        for key in keys(testdisaggregatetranslateplan.translatingnodes[i])
            @test sort(testdisaggregatetranslateplan.translatingnodes[i][key]) ==
                sort(ptestdisaggregatetranslateplan.translatingnodes[i][key])
        end
    end

    @test ptestaggregatetranslateplan.levels == testaggregatetranslateplan.levels
    @test ptestaggregatetranslateplan.rootoffset == testaggregatetranslateplan.rootoffset
    @test ptestaggregatetranslateplan.tree == testaggregatetranslateplan.tree
    for level in testaggregatetranslateplan.levels
        @test sort(H2Trees.aggregationnodes(ptestaggregatetranslateplan, level)) ==
            sort(H2Trees.aggregationnodes(testaggregatetranslateplan, level))
    end
    for i in eachindex(testaggregatetranslateplan.receivingnodes)
        for key in keys(testaggregatetranslateplan.receivingnodes[i])
            @test sort(testaggregatetranslateplan.receivingnodes[i][key]) ==
                sort(ptestaggregatetranslateplan.receivingnodes[i][key])
        end
    end

    @test ptrialdisaggregateplan.levels == trialdisaggregateplan.levels
    @test ptrialdisaggregateplan.storenode == trialdisaggregateplan.storenode
    @test ptrialdisaggregateplan.tree == trialdisaggregateplan.tree
    @test ptrialdisaggregateplan.rootoffset == trialdisaggregateplan.rootoffset
    for level in trialdisaggregateplan.levels
        @test sort(H2Trees.disaggregationnodes(ptrialdisaggregateplan, level)) ==
            sort(H2Trees.disaggregationnodes(trialdisaggregateplan, level))
    end

    trialaggregatetranslateplan = H2Trees.AggregateTranslatePlan(tree, TFIterator(tree))
    testdisaggregateplan = H2Trees.DisaggregatePlan(tree, aggregatenode(tree))
    testaggregateplan = H2Trees.AggregatePlan(tree, aggregatenode(tree))
    trialdisaggregatetranslateplan = H2Trees.DisaggregateTranslatePlan(
        tree, TFIterator(tree)
    )
    plans = H2Trees._buildgalerkinplans(
        tree, aggregatenode, TFIterator, H2Trees.AggregateTranslateMode()
    )

    ptrialaggregatetranslateplan = plans.trialaggregationplan
    ptestdisaggregateplan = plans.testdisaggregationplan
    ptestaggregateplan = plans.testaggregationplan
    ptrialdisaggregatetranslateplan = plans.trialdisaggregationplan

    @test ptrialaggregatetranslateplan.levels == trialaggregatetranslateplan.levels
    @test ptrialaggregatetranslateplan.rootoffset == trialaggregatetranslateplan.rootoffset
    @test ptrialaggregatetranslateplan.tree == trialaggregatetranslateplan.tree
    for level in trialaggregatetranslateplan.levels
        @test sort(H2Trees.aggregationnodes(ptrialaggregatetranslateplan, level)) ==
            sort(H2Trees.aggregationnodes(trialaggregatetranslateplan, level))
    end
    for i in eachindex(trialaggregatetranslateplan.receivingnodes)
        for key in keys(trialaggregatetranslateplan.receivingnodes[i])
            @test sort(trialaggregatetranslateplan.receivingnodes[i][key]) ==
                sort(ptrialaggregatetranslateplan.receivingnodes[i][key])
        end
    end

    @test ptestdisaggregateplan.levels == testdisaggregateplan.levels
    @test ptestdisaggregateplan.storenode == testdisaggregateplan.storenode
    @test ptestdisaggregateplan.tree == testdisaggregateplan.tree
    @test ptestdisaggregateplan.rootoffset == testdisaggregateplan.rootoffset
    for level in testdisaggregateplan.levels
        @test sort(H2Trees.disaggregationnodes(ptestdisaggregateplan, level)) ==
            sort(H2Trees.disaggregationnodes(testdisaggregateplan, level))
    end

    @test ptestaggregateplan.levels == testaggregateplan.levels
    @test ptestaggregateplan.storenode == testaggregateplan.storenode
    @test ptestaggregateplan.tree == testaggregateplan.tree
    @test ptestaggregateplan.rootoffset == testaggregateplan.rootoffset
    for level in testaggregateplan.levels
        @test sort(H2Trees.aggregationnodes(ptestaggregateplan, level)) ==
            sort(H2Trees.aggregationnodes(testaggregateplan, level))
    end

    @test ptrialdisaggregatetranslateplan.levels == trialdisaggregatetranslateplan.levels
    @test ptrialdisaggregatetranslateplan.isdisaggregationnode ==
        trialdisaggregatetranslateplan.isdisaggregationnode
    @test ptrialdisaggregatetranslateplan.tree == trialdisaggregatetranslateplan.tree
    @test ptrialdisaggregatetranslateplan.rootoffset ==
        trialdisaggregatetranslateplan.rootoffset
    for level in trialdisaggregatetranslateplan.levels
        @test sort(H2Trees.disaggregationnodes(ptrialdisaggregatetranslateplan, level)) ==
            sort(H2Trees.disaggregationnodes(trialdisaggregatetranslateplan, level))
    end
    @test sort(collect(keys(trialdisaggregatetranslateplan.translatingnodes))) ==
        sort(collect(keys(ptrialdisaggregatetranslateplan.translatingnodes)))
    for i in eachindex(trialdisaggregatetranslateplan.translatingnodes)
        for key in keys(trialdisaggregatetranslateplan.translatingnodes[i])
            @test sort(trialdisaggregatetranslateplan.translatingnodes[i][key]) ==
                sort(ptrialdisaggregatetranslateplan.translatingnodes[i][key])
        end
    end
end

@testset "Petrov Plan" begin
    λ = 1.0

    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere6.in")
    )
    m1 = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter14.in")
    )
    m2 = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter15.in")
    )
    m3 = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter16.in")
    )

    for mx in [m, m2, m3]
        X = raviartthomas(mx)
        for my in [m, m2, m3]
            Y = raviartthomas(my)
            tree = buildtree(
                X,
                Y;
                builder=BlockTreeBuilder(;
                    test=TwoNTreeBuilder(; minhalfsize=λ / 10, minvalues=0),
                    trial=TwoNTreeBuilder(; minhalfsize=λ / 10, minvalues=0),
                ),
            )
            TFIterator = H2Trees.TranslatingNodesIterator(; isnear=H2Trees.isnear())
            aggregatenode = H2Trees.istranslatingnode(; TranslatingNodesIterator=TFIterator)

            testtree = H2Trees.testtree(tree)
            trialtree = H2Trees.trialtree(tree)

            trialaggregateplan = H2Trees.AggregatePlan(
                trialtree,
                H2Trees.PetrovAggregationFunctor(aggregatenode, tree, testtree, trialtree),
            )
            testdisaggregatetranslateplan = H2Trees.DisaggregateTranslatePlan(
                testtree, trialtree, H2Trees.TranslatingNodesIterator
            )
            testaggregatetranslateplan = H2Trees.AggregateTranslatePlan(
                testtree, H2Trees.PetrovDisaggregationFunctor(TFIterator, tree, trialtree)
            )
            trialdisaggregateplan = H2Trees.DisaggregatePlan(
                trialtree,
                H2Trees.PetrovAggregationFunctor(aggregatenode, tree, testtree, trialtree),
            )

            aggregateallplans = H2Trees._buildpetrovplans(
                tree, H2Trees.aggregateallnodes(), TFIterator, H2Trees.AggregateMode()
            )
            @test sum(aggregateallplans.trialaggregationplan.storenode) ==
                length(trialtree.nodes)

            aggregaterootplans = H2Trees._buildpetrovplans(
                tree, H2Trees.aggregaterootonly(), TFIterator, H2Trees.AggregateMode()
            )
            @test H2Trees.aggregaterootonly()(tree)(
                testtree, trialtree, H2Trees.root(trialtree)
            )
            @test sum(aggregaterootplans.trialaggregationplan.storenode) == 1
            aggregaterootplans.trialaggregationplan.storenode[H2Trees.root(trialtree)] == 1
            plans = H2Trees._buildpetrovplans(
                tree, aggregatenode, TFIterator, H2Trees.AggregateMode()
            )

            @test plans isa H2Trees.PlanSet
            @test H2Trees.ispetrov(plans)
            @test !H2Trees.isgalerkin(plans)
            @test H2Trees.tree(plans) === tree
            @test H2Trees.testtree(plans) === testtree
            @test H2Trees.trialtree(plans) === trialtree
            @test H2Trees.mintranslationlevel(plans) == plans.mintranslationlevel
            @test hasproperty(plans, :tree)
            @test hasproperty(plans, :family)
            @test H2Trees.planfamily(plans) isa H2Trees.PetrovPlanFamily
            @test :mintranslationlevel in keys(plans)
            @test plans[:mintranslationlevel] == plans.mintranslationlevel
            @test NamedTuple(plans) == (
                testaggregationplan=plans.testaggregationplan,
                trialaggregationplan=plans.trialaggregationplan,
                testdisaggregationplan=plans.testdisaggregationplan,
                trialdisaggregationplan=plans.trialdisaggregationplan,
                relevantlevels=plans.relevantlevels,
                mintranslationlevel=plans.mintranslationlevel,
            )
            @test (; plans...) == NamedTuple(plans)
            @test merge((before=:before,), plans) ==
                merge((before=:before,), NamedTuple(plans))
            @test merge(plans, (after=:after,)) == merge(NamedTuple(plans), (after=:after,))
            @test H2Trees.ownedtree(plans.testdisaggregationplan) === testtree
            @test H2Trees.receivingtree(plans, plans.testdisaggregationplan) === testtree
            @test H2Trees.translatingtree(plans, plans.testdisaggregationplan) === trialtree
            @test H2Trees.receivingtree(tree, plans.testaggregationplan) === trialtree
            @test H2Trees.translatingtree(tree, plans.testaggregationplan) === testtree

            builtplans = H2Trees.buildplans(
                tree;
                builder=H2Trees.PlanBuilder(;
                    aggregationmode=H2Trees.AggregateMode(), isnear=H2Trees.isnear()
                ),
            )
            @test builtplans isa H2Trees.PlanSet
            @test H2Trees.ispetrov(builtplans)
            @test H2Trees.relevantlevels(builtplans) == H2Trees.relevantlevels(plans)

            ptrialaggregateplan = plans.trialaggregationplan
            ptestdisaggregatetranslateplan = plans.testdisaggregationplan
            ptestaggregatetranslateplan = plans.testaggregationplan
            ptrialdisaggregateplan = plans.trialdisaggregationplan

            @test ptrialaggregateplan.levels == trialaggregateplan.levels
            @test ptrialaggregateplan.storenode == trialaggregateplan.storenode
            @test ptrialaggregateplan.tree == trialaggregateplan.tree
            @test ptrialaggregateplan.rootoffset == trialaggregateplan.rootoffset
            for level in trialaggregateplan.levels
                @test sort(H2Trees.aggregationnodes(ptrialaggregateplan, level)) ==
                    sort(H2Trees.aggregationnodes(trialaggregateplan, level))
            end

            @test ptestdisaggregatetranslateplan.levels ==
                testdisaggregatetranslateplan.levels
            @test ptestdisaggregatetranslateplan.isdisaggregationnode ==
                testdisaggregatetranslateplan.isdisaggregationnode
            @test ptestdisaggregatetranslateplan.tree == testdisaggregatetranslateplan.tree
            @test ptestdisaggregatetranslateplan.rootoffset ==
                testdisaggregatetranslateplan.rootoffset
            for level in testdisaggregatetranslateplan.levels
                @test sort(
                    H2Trees.disaggregationnodes(ptestdisaggregatetranslateplan, level)
                ) == sort(
                    H2Trees.disaggregationnodes(testdisaggregatetranslateplan, level)
                )
            end
            @test sort(collect(keys(testdisaggregatetranslateplan.translatingnodes))) ==
                sort(collect(keys(ptestdisaggregatetranslateplan.translatingnodes)))
            for i in eachindex(testdisaggregatetranslateplan.translatingnodes)
                for key in keys(testdisaggregatetranslateplan.translatingnodes[i])
                    @test sort(testdisaggregatetranslateplan.translatingnodes[i][key]) ==
                        sort(ptestdisaggregatetranslateplan.translatingnodes[i][key])
                end
            end

            @test ptestaggregatetranslateplan.levels == testaggregatetranslateplan.levels
            @test ptestaggregatetranslateplan.rootoffset ==
                testaggregatetranslateplan.rootoffset
            @test ptestaggregatetranslateplan.tree == testaggregatetranslateplan.tree
            for level in testaggregatetranslateplan.levels
                @test sort(H2Trees.aggregationnodes(ptestaggregatetranslateplan, level)) ==
                    sort(H2Trees.aggregationnodes(testaggregatetranslateplan, level))
            end
            for i in eachindex(testaggregatetranslateplan.receivingnodes)
                for key in keys(testaggregatetranslateplan.receivingnodes[i])
                    @test sort(testaggregatetranslateplan.receivingnodes[i][key]) ==
                        sort(ptestaggregatetranslateplan.receivingnodes[i][key])
                end
            end

            @test ptrialdisaggregateplan.levels == trialdisaggregateplan.levels
            @test ptrialdisaggregateplan.storenode == trialdisaggregateplan.storenode
            @test ptrialdisaggregateplan.tree == trialdisaggregateplan.tree
            @test ptrialdisaggregateplan.rootoffset == trialdisaggregateplan.rootoffset
            for level in trialdisaggregateplan.levels
                @test sort(H2Trees.disaggregationnodes(ptrialdisaggregateplan, level)) ==
                    sort(H2Trees.disaggregationnodes(trialdisaggregateplan, level))
            end

            trialaggregatetranslateplan = H2Trees.AggregateTranslatePlan(
                trialtree, H2Trees.PetrovDisaggregationFunctor(TFIterator, tree, testtree)
            )
            testdisaggregateplan = H2Trees.DisaggregatePlan(
                testtree,
                H2Trees.PetrovAggregationFunctor(aggregatenode, tree, trialtree, testtree),
            )
            testaggregateplan = H2Trees.AggregatePlan(
                testtree,
                H2Trees.PetrovAggregationFunctor(aggregatenode, tree, trialtree, testtree),
            )
            trialdisaggregatetranslateplan = H2Trees.DisaggregateTranslatePlan(
                trialtree, testtree, H2Trees.TranslatingNodesIterator
            )

            plans = H2Trees._buildpetrovplans(
                tree, aggregatenode, TFIterator, H2Trees.AggregateTranslateMode()
            )
            ptrialaggregatetranslateplan = plans.trialaggregationplan
            ptestdisaggregateplan = plans.testdisaggregationplan
            ptestaggregateplan = plans.testaggregationplan
            ptrialdisaggregatetranslateplan = plans.trialdisaggregationplan

            @test H2Trees.receivingtree(plans, plans.trialaggregationplan) === testtree
            @test H2Trees.translatingtree(plans, plans.trialaggregationplan) === trialtree
            @test H2Trees.receivingtree(tree, plans.trialdisaggregationplan) === trialtree
            @test H2Trees.translatingtree(tree, plans.trialdisaggregationplan) === testtree

            @test ptrialaggregatetranslateplan.levels == trialaggregatetranslateplan.levels
            @test ptrialaggregatetranslateplan.rootoffset ==
                trialaggregatetranslateplan.rootoffset
            @test ptrialaggregatetranslateplan.tree == trialaggregatetranslateplan.tree
            for level in trialaggregatetranslateplan.levels
                @test sort(H2Trees.aggregationnodes(ptrialaggregatetranslateplan, level)) ==
                    sort(H2Trees.aggregationnodes(trialaggregatetranslateplan, level))
            end
            for i in eachindex(trialaggregatetranslateplan.receivingnodes)
                for key in keys(trialaggregatetranslateplan.receivingnodes[i])
                    @test sort(trialaggregatetranslateplan.receivingnodes[i][key]) ==
                        sort(ptrialaggregatetranslateplan.receivingnodes[i][key])
                end
            end

            @test ptestdisaggregateplan.levels == testdisaggregateplan.levels
            @test ptestdisaggregateplan.storenode == testdisaggregateplan.storenode
            @test ptestdisaggregateplan.tree == testdisaggregateplan.tree
            @test ptestdisaggregateplan.rootoffset == testdisaggregateplan.rootoffset
            for level in testdisaggregateplan.levels
                @test sort(H2Trees.disaggregationnodes(ptestdisaggregateplan, level)) ==
                    sort(H2Trees.disaggregationnodes(testdisaggregateplan, level))
            end

            @test ptestaggregateplan.levels == testaggregateplan.levels
            @test ptestaggregateplan.storenode == testaggregateplan.storenode
            @test ptestaggregateplan.tree == testaggregateplan.tree
            @test ptestaggregateplan.rootoffset == testaggregateplan.rootoffset
            for level in testaggregateplan.levels
                @test sort(H2Trees.aggregationnodes(ptestaggregateplan, level)) ==
                    sort(H2Trees.aggregationnodes(testaggregateplan, level))
            end

            @test ptrialdisaggregatetranslateplan.levels ==
                trialdisaggregatetranslateplan.levels
            @test ptrialdisaggregatetranslateplan.isdisaggregationnode ==
                trialdisaggregatetranslateplan.isdisaggregationnode
            @test ptrialdisaggregatetranslateplan.tree ==
                trialdisaggregatetranslateplan.tree
            @test ptrialdisaggregatetranslateplan.rootoffset ==
                trialdisaggregatetranslateplan.rootoffset
            for level in trialdisaggregatetranslateplan.levels
                @test sort(
                    H2Trees.disaggregationnodes(ptrialdisaggregatetranslateplan, level)
                ) == sort(
                    H2Trees.disaggregationnodes(trialdisaggregatetranslateplan, level)
                )
            end
            @test sort(collect(keys(trialdisaggregatetranslateplan.translatingnodes))) ==
                sort(collect(keys(ptrialdisaggregatetranslateplan.translatingnodes)))
            for i in eachindex(trialdisaggregatetranslateplan.translatingnodes)
                for key in keys(trialdisaggregatetranslateplan.translatingnodes[i])
                    @test sort(trialdisaggregatetranslateplan.translatingnodes[i][key]) ==
                        sort(ptrialdisaggregatetranslateplan.translatingnodes[i][key])
                end
            end
        end
    end
end
