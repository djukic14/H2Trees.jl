using Test
using BEAST, CompScienceMeshes
using StaticArrays
using LinearAlgebra
using H2Trees

@testset "Galerkin Plan" begin
    λ = 1.0
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere6.in")
    )
    X = raviartthomas(m)
    tree = TwoNTree(X, λ / 10)

    TFIterator = H2Trees.TranslatingNodesIterator(; isnear=H2Trees.isnear())
    aggregatenode = H2Trees.istranslatingnode(; TranslatingNodesIterator=TFIterator)

    trialaggregateplan = H2Trees.AggregatePlan(tree, aggregatenode(tree))
    testdisaggregatetranslateplan = H2Trees.DisaggregateTranslatePlan(
        tree, TFIterator(tree)
    )
    testaggregatetranslateplan = H2Trees.AggregateTranslatePlan(tree, TFIterator(tree))
    trialdisaggregateplan = H2Trees.DisaggregatePlan(tree, aggregatenode(tree))

    plans = H2Trees.galerkinplans(tree, aggregatenode, TFIterator, H2Trees.AggregateMode())

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
    plans = H2Trees.galerkinplans(
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
            tree = TwoNTree(X, Y, λ / 10)
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
            plans = H2Trees.petrovplans(
                tree, aggregatenode, TFIterator, H2Trees.AggregateMode()
            )

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

            plans = H2Trees.petrovplans(
                tree, aggregatenode, TFIterator, H2Trees.AggregateTranslateMode()
            )
            ptrialaggregatetranslateplan = plans.trialaggregationplan
            ptestdisaggregateplan = plans.testdisaggregationplan
            ptestaggregateplan = plans.testaggregationplan
            ptrialdisaggregatetranslateplan = plans.trialdisaggregationplan

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
