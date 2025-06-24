using Test
using BEAST, CompScienceMeshes
using StaticArrays
using LinearAlgebra
using H2Trees

@testset "Adjoint Plans" begin
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

    for mx in [m, m1, m2, m3]
        X = raviartthomas(mx)
        for my in [m, m1, m2, m3]
            Y = raviartthomas(my)
            tree = TwoNTree(X, Y, λ / 10)
            TFIterator = H2Trees.TranslatingNodesIterator(; isnear=H2Trees.isnear())
            aggregatenode = H2Trees.istranslatingnode(; TranslatingNodesIterator=TFIterator)

            testtree = H2Trees.testtree(tree)
            trialtree = H2Trees.trialtree(tree)

            @test_throws ErrorException H2Trees.DisaggregateTranslatePlan(
                tree, H2Trees.TranslatingNodesIterator
            )

            @test_throws ErrorException H2Trees.AggregateTranslatePlan(
                tree, H2Trees.TranslatingNodesIterator
            )

            trialaggregateplan = H2Trees.AggregatePlan(
                trialtree,
                H2Trees.PetrovAggregationFunctor(aggregatenode, tree, testtree, trialtree),
            )
            testdisaggregatetranslateplan = H2Trees.DisaggregateTranslatePlan(
                testtree, trialtree, H2Trees.TranslatingNodesIterator
            )
            @test testdisaggregatetranslateplan[1, 1] == Int[]

            atestaggregatetranslateplan, atrialdisaggregateplan = H2Trees.adjointplans(
                trialaggregateplan, testdisaggregatetranslateplan
            )

            testaggregatetranslateplan = H2Trees.AggregateTranslatePlan(
                testtree, H2Trees.PetrovDisaggregationFunctor(TFIterator, tree, trialtree)
            )

            trialdisaggregateplan = H2Trees.DisaggregatePlan(
                trialtree,
                H2Trees.PetrovAggregationFunctor(aggregatenode, tree, testtree, trialtree),
            )

            @test H2Trees.translatingplan(
                trialaggregateplan, testdisaggregatetranslateplan
            ) == testdisaggregatetranslateplan
            @test H2Trees.translatingplan(
                testaggregatetranslateplan, trialdisaggregateplan
            ) == testaggregatetranslateplan
            @test_throws ErrorException H2Trees.translatingplan(
                trialaggregateplan, trialdisaggregateplan
            )
            @test_throws ErrorException H2Trees.translatingplan(
                testaggregatetranslateplan, testdisaggregatetranslateplan
            )

            @test H2Trees.aggregationlevels(atestaggregatetranslateplan) ==
                H2Trees.aggregationlevels(testaggregatetranslateplan)
            @test H2Trees.rootoffset(atestaggregatetranslateplan) ==
                H2Trees.rootoffset(testaggregatetranslateplan)
            @test H2Trees.tree(atestaggregatetranslateplan) ==
                H2Trees.tree(testaggregatetranslateplan)

            for level in H2Trees.aggregationlevels(testaggregatetranslateplan)
                @test sort(H2Trees.aggregationnodes(atestaggregatetranslateplan, level)) ==
                    sort(H2Trees.aggregationnodes(testaggregatetranslateplan, level))
            end

            for i in eachindex(H2Trees.receivingnodes(testaggregatetranslateplan))
                rnodes = H2Trees.receivingnodes(testaggregatetranslateplan)[i]
                arnodes = H2Trees.receivingnodes(atestaggregatetranslateplan)[i]

                for k in keys(rnodes)
                    @test sort(rnodes[k]) == sort(arnodes[k])

                    @test sort(
                        H2Trees.translatingnodes(
                            testaggregatetranslateplan, k, H2Trees.level(trialtree, k)
                        ),
                    ) == sort(
                        H2Trees.translatingnodes(
                            atestaggregatetranslateplan, k, H2Trees.level(trialtree, k)
                        ),
                    )
                end
            end

            @test testaggregatetranslateplan[100000, 3] == Int[]

            @test H2Trees.storenode(trialdisaggregateplan) ==
                H2Trees.storenode(atrialdisaggregateplan)
            @test H2Trees.disaggregationlevels(trialdisaggregateplan) ==
                H2Trees.disaggregationlevels(atrialdisaggregateplan)
            @test H2Trees.rootoffset(trialdisaggregateplan) ==
                H2Trees.rootoffset(atrialdisaggregateplan)
            @test H2Trees.tree(trialdisaggregateplan) ==
                H2Trees.tree(atrialdisaggregateplan)

            for level in H2Trees.disaggregationlevels(trialdisaggregateplan)
                @test sort(H2Trees.disaggregationnodes(trialdisaggregateplan, level)) ==
                    sort(H2Trees.disaggregationnodes(atrialdisaggregateplan, level))
            end

            trialaggregatetranslateplan = H2Trees.AggregateTranslatePlan(
                trialtree, H2Trees.PetrovDisaggregationFunctor(TFIterator, tree, testtree)
            )
            testdisaggregateplan = H2Trees.DisaggregatePlan(
                testtree,
                H2Trees.PetrovAggregationFunctor(aggregatenode, tree, trialtree, testtree),
            )
            atestaggregateplan, atrialdisaggregatetranslateplan = H2Trees.adjointplans(
                trialaggregatetranslateplan, testdisaggregateplan
            )

            testaggregateplan = H2Trees.AggregatePlan(
                testtree,
                H2Trees.PetrovAggregationFunctor(aggregatenode, tree, trialtree, testtree),
            )
            trialdisaggregatetranslateplan = H2Trees.DisaggregateTranslatePlan(
                trialtree, testtree, H2Trees.TranslatingNodesIterator
            )

            @test H2Trees.aggregationlevels(atestaggregateplan) ==
                H2Trees.aggregationlevels(testaggregateplan)
            @test H2Trees.rootoffset(atestaggregateplan) ==
                H2Trees.rootoffset(testaggregateplan)
            @test H2Trees.tree(atestaggregateplan) == H2Trees.tree(testaggregateplan)
            @test H2Trees.storenode(atestaggregateplan) ==
                H2Trees.storenode(testaggregateplan)

            for level in H2Trees.aggregationlevels(testaggregateplan)
                @test sort(H2Trees.aggregationnodes(atestaggregateplan, level)) ==
                    sort(H2Trees.aggregationnodes(testaggregateplan, level))
            end

            @test H2Trees.disaggregationlevels(atrialdisaggregatetranslateplan) ==
                H2Trees.disaggregationlevels(trialdisaggregatetranslateplan)
            @test H2Trees.rootoffset(atrialdisaggregatetranslateplan) ==
                H2Trees.rootoffset(trialdisaggregatetranslateplan)
            @test H2Trees.tree(atrialdisaggregatetranslateplan) ==
                H2Trees.tree(trialdisaggregatetranslateplan)
            @test atrialdisaggregatetranslateplan.isdisaggregationnode ==
                trialdisaggregatetranslateplan.isdisaggregationnode

            for level in H2Trees.disaggregationlevels(trialdisaggregatetranslateplan)
                @test sort(
                    H2Trees.disaggregationnodes(atrialdisaggregatetranslateplan, level)
                ) == sort(
                    H2Trees.disaggregationnodes(trialdisaggregatetranslateplan, level)
                )
            end

            @test trialaggregatetranslateplan[100000, 3] == Int[]
            rnode = first(keys(trialdisaggregatetranslateplan.translatingnodes[1]))
            @test trialdisaggregatetranslateplan[rnode] ==
                trialdisaggregatetranslateplan.translatingnodes[1][rnode] ==
                H2Trees.translatingnodes(trialdisaggregatetranslateplan, rnode)

            for i in eachindex(H2Trees.translatingnodes(trialdisaggregatetranslateplan))
                tfnodes = H2Trees.translatingnodes(trialdisaggregatetranslateplan)[i]
                atfnodes = H2Trees.translatingnodes(atrialdisaggregatetranslateplan)[i]

                for key in keys(tfnodes)
                    @test sort(tfnodes[key]) == sort(atfnodes[key])
                end
            end
        end
    end
end
