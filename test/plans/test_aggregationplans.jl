using Test
using BEAST, CompScienceMeshes
using StaticArrays
using LinearAlgebra
using H2Trees

@testset "Aggregation Plan" begin
    λ = 1.0
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere6.in")
    )
    X = raviartthomas(m)
    root = 2
    tree = TwoNTree(X, λ / 10; root=root)

    aggregationplan = H2Trees.AggregatePlan(tree, node -> true)
    for node in H2Trees.DepthFirstIterator(tree)
        @test H2Trees.storenode(aggregationplan, node)
    end

    @test aggregationplan.levels == reverse(H2Trees.levels(tree))
    aggregationnodes = reverse(H2Trees.nodesatlevel(tree))
    for (i, nodes) in enumerate(aggregationnodes)
        @test sort(aggregationplan.nodes[i]) == sort(nodes)
    end

    aggregationplan = H2Trees.AggregatePlan(tree, node -> ifelse(node == root, true, false))
    @test aggregationplan.levels == reverse(H2Trees.levels(tree))
    for (i, nodes) in enumerate(aggregationnodes)
        @test sort(aggregationplan.nodes[i]) == sort(nodes)
    end

    for node in H2Trees.DepthFirstIterator(tree)
        if node == root
            @test H2Trees.storenode(aggregationplan, node)
        else
            @test !H2Trees.storenode(aggregationplan, node)
        end
    end

    for level in H2Trees.aggregationlevels(aggregationplan)
        @test unique(H2Trees.aggregationnodes(aggregationplan, level)) ==
            H2Trees.aggregationnodes(aggregationplan, level)
        for node in H2Trees.aggregationnodes(aggregationplan, level)
            @test H2Trees.level(tree, node) == level

            aggregationnode = false

            if H2Trees.root(tree) == node
                @test H2Trees.storenode(aggregationplan, node)
                aggregationnode = true
            else
                for parent in H2Trees.ParentUpwardsIterator(tree, node)
                    if H2Trees.root(tree) == parent
                        aggregationnode = true
                        break
                    end
                end
            end

            @test aggregationnode
        end
    end

    aggregationplan = H2Trees.AggregatePlan(
        tree, node -> H2Trees.istranslatingnode(tree, node)
    )

    for node in H2Trees.DepthFirstIterator(tree, root)
        if H2Trees.istranslatingnode(tree, node)
            @test H2Trees.storenode(aggregationplan, node)
        else
            @test !H2Trees.storenode(aggregationplan, node)
        end
    end

    for level in H2Trees.aggregationlevels(aggregationplan)
        @test unique(H2Trees.aggregationnodes(aggregationplan, level)) ==
            H2Trees.aggregationnodes(aggregationplan, level)
        for node in H2Trees.aggregationnodes(aggregationplan, level)
            @test H2Trees.level(tree, node) == level

            aggregationnode = false

            if H2Trees.istranslatingnode(tree, node)
                @test H2Trees.storenode(aggregationplan, node)
                aggregationnode = true
            else
                for parent in H2Trees.ParentUpwardsIterator(tree, node)
                    if H2Trees.istranslatingnode(tree, parent)
                        aggregationnode = true
                        break
                    end
                end
            end

            @test aggregationnode
        end
    end

    @test H2Trees.rootoffset(aggregationplan) == 1

    @test H2Trees.mintranslationlevel(tree) ==
        H2Trees.minaggregationlevel(aggregationplan) ==
        3

    m1 = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter8.in")
    )
    m2 = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter9.in")
    )
    m3 = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter10.in")
    )

    for mx in [m, m1, m2, m3]
        X = raviartthomas(mx)
        for my in [m, m1, m2, m3]
            Y = raviartthomas(my)

            tree = TwoNTree(X, Y, λ / 10)
            testtree = H2Trees.testtree(tree)
            trialtree = H2Trees.trialtree(tree)

            @test_throws ErrorException H2Trees.AggregatePlan(tree, node -> true)

            trialaggregationplan = H2Trees.AggregatePlan(
                trialtree, node -> H2Trees.istranslatingnode(testtree, trialtree, node)
            )

            for trialnode in H2Trees.DepthFirstIterator(trialtree, H2Trees.root(trialtree))
                if H2Trees.istranslatingnode(testtree, trialtree, trialnode)
                    @test H2Trees.storenode(trialaggregationplan, trialnode)
                else
                    @test !H2Trees.storenode(trialaggregationplan, trialnode)
                end
            end

            for level in H2Trees.aggregationlevels(trialaggregationplan)
                @test unique(H2Trees.aggregationnodes(trialaggregationplan, level)) ==
                    H2Trees.aggregationnodes(trialaggregationplan, level)

                for trialnode in H2Trees.aggregationnodes(trialaggregationplan, level)
                    @test H2Trees.level(trialtree, trialnode) == level

                    aggregationnode = false
                    if H2Trees.istranslatingnode(testtree, trialtree, trialnode)
                        @test H2Trees.storenode(trialaggregationplan, trialnode)
                        aggregationnode = true
                    else
                        for parent in H2Trees.ParentUpwardsIterator(trialtree, trialnode)
                            if H2Trees.istranslatingnode(testtree, trialtree, parent)
                                aggregationnode = true
                                break
                            end
                        end
                    end

                    @test aggregationnode
                end
            end

            @test H2Trees.rootoffset(trialaggregationplan) == H2Trees.root(trialtree) - 1

            testaggregationplan = H2Trees.AggregatePlan(
                testtree, node -> H2Trees.istranslatingnode(trialtree, testtree, node)
            )

            for testnode in H2Trees.DepthFirstIterator(testtree, H2Trees.root(testtree))
                if H2Trees.istranslatingnode(trialtree, testtree, testnode)
                    @test H2Trees.storenode(testaggregationplan, testnode)
                else
                    @test !H2Trees.storenode(testaggregationplan, testnode)
                end
            end

            for level in H2Trees.aggregationlevels(testaggregationplan)
                @test unique(H2Trees.aggregationnodes(testaggregationplan, level)) ==
                    H2Trees.aggregationnodes(testaggregationplan, level)

                for testnode in H2Trees.aggregationnodes(testaggregationplan, level)
                    @test H2Trees.level(testtree, testnode) == level

                    aggregationnode = false
                    if H2Trees.istranslatingnode(trialtree, testtree, testnode)
                        @test H2Trees.storenode(testaggregationplan, testnode)
                        aggregationnode = true
                    else
                        for parent in H2Trees.ParentUpwardsIterator(testtree, testnode)
                            if H2Trees.istranslatingnode(trialtree, testtree, parent)
                                aggregationnode = true
                                break
                            end
                        end
                    end

                    @test aggregationnode
                end
            end

            @test H2Trees.rootoffset(testaggregationplan) == H2Trees.root(testtree) - 1
        end
    end
end

@testset "Aggregate Translate Plan" begin
    λ = 1.0
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere6.in")
    )
    X = raviartthomas(m)
    tree = TwoNTree(X, λ / 10)

    TFIterator = H2Trees.TranslatingNodesIterator(; isnear=H2Trees.isnear())(tree)
    aggregatetranslateplan = H2Trees.AggregateTranslatePlan(tree, TFIterator)
    disaggregationplan = H2Trees.DisaggregateTranslatePlan(tree, TFIterator)

    aggregateplan = H2Trees.AggregatePlan(
        tree, node -> H2Trees.istranslatingnode(tree, node)
    )

    @test H2Trees.aggregationlevels(aggregatetranslateplan) ==
        H2Trees.aggregationlevels(aggregateplan)

    @test H2Trees.rootoffset(aggregatetranslateplan) == H2Trees.rootoffset(aggregateplan)
    @test H2Trees.minaggregationlevel(aggregatetranslateplan) ==
        H2Trees.minaggregationlevel(aggregateplan)
    @test H2Trees.mintranslationlevel(aggregatetranslateplan) ==
        H2Trees.mintranslationlevel(disaggregationplan)

    for level in H2Trees.aggregationlevels(aggregatetranslateplan)
        @test sort(H2Trees.aggregationnodes(aggregatetranslateplan, level)) ==
            sort(H2Trees.aggregationnodes(aggregateplan, level))
    end

    N = length(H2Trees.aggregationlevels(aggregatetranslateplan))
    for i in eachindex(disaggregationplan.nodes)
        tfd = H2Trees.translatingnodes(disaggregationplan)[i]
        tfa = H2Trees.receivingnodes(aggregatetranslateplan)[N - i + 1]

        for k in keys(tfd)
            @test sort(tfd[k]) == sort(tfa[k])
        end
    end

    for level in H2Trees.aggregationlevels(aggregatetranslateplan)
        @test sort(H2Trees.aggregationnodes(aggregatetranslateplan, level)) ==
            sort(H2Trees.disaggregationnodes(disaggregationplan, level))
    end

    N = length(H2Trees.aggregationlevels(aggregatetranslateplan))
    for i in eachindex(disaggregationplan.nodes)
        tfd = H2Trees.translatingnodes(disaggregationplan)[i]
        tfa = H2Trees.receivingnodes(aggregatetranslateplan)[N - i + 1]

        for k in keys(tfd)
            @test sort(tfd[k]) == sort(tfa[k])
        end
    end

    m1 = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter8.in")
    )
    m2 = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter9.in")
    )
    m3 = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter10.in")
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

            @test_throws ErrorException H2Trees.AggregateTranslatePlan(
                tree, H2Trees.TranslatingNodesIterator
            )

            trialaggregatetranslateplan = H2Trees.AggregateTranslatePlan(
                trialtree, H2Trees.PetrovDisaggregationFunctor(TFIterator, tree, testtree)
            )

            testaggregatetranslateplan = H2Trees.AggregateTranslatePlan(
                testtree, H2Trees.PetrovDisaggregationFunctor(TFIterator, tree, trialtree)
            )

            trialaggregateplan = H2Trees.AggregatePlan(
                trialtree,
                H2Trees.PetrovAggregationFunctor(aggregatenode, tree, testtree, trialtree),
            )

            testaggregateplan = H2Trees.AggregatePlan(
                testtree,
                H2Trees.PetrovAggregationFunctor(aggregatenode, tree, trialtree, testtree),
            )

            trialdisaggregationplan = H2Trees.DisaggregateTranslatePlan(
                trialtree, testtree, H2Trees.TranslatingNodesIterator
            )
            testdisaggregationplan = H2Trees.DisaggregateTranslatePlan(
                testtree, trialtree, H2Trees.TranslatingNodesIterator
            )

            @test_throws ErrorException H2Trees.DisaggregateTranslatePlan(
                tree, H2Trees.TranslatingNodesIterator
            )

            @test H2Trees.aggregationlevels(trialaggregatetranslateplan) ==
                H2Trees.aggregationlevels(trialaggregateplan)
            @test H2Trees.aggregationlevels(testaggregatetranslateplan) ==
                H2Trees.aggregationlevels(testaggregateplan)

            @test H2Trees.rootoffset(trialaggregatetranslateplan) ==
                H2Trees.rootoffset(trialaggregateplan)
            @test H2Trees.rootoffset(testaggregatetranslateplan) ==
                H2Trees.rootoffset(testaggregateplan)

            @test H2Trees.minaggregationlevel(trialaggregatetranslateplan) ==
                H2Trees.minaggregationlevel(trialaggregateplan)

            @test H2Trees.minaggregationlevel(testaggregatetranslateplan) ==
                H2Trees.minaggregationlevel(testaggregateplan)

            @test H2Trees.mintranslationlevel(trialaggregatetranslateplan) ==
                H2Trees.mintranslationlevel(trialdisaggregationplan)
            @test H2Trees.mintranslationlevel(testaggregatetranslateplan) ==
                H2Trees.mintranslationlevel(testdisaggregationplan)

            for level in H2Trees.aggregationlevels(trialaggregatetranslateplan)
                @test sort(H2Trees.aggregationnodes(trialaggregatetranslateplan, level)) ==
                    sort(H2Trees.aggregationnodes(trialaggregateplan, level))
            end

            for level in H2Trees.aggregationlevels(testaggregatetranslateplan)
                @test sort(H2Trees.aggregationnodes(testaggregatetranslateplan, level)) ==
                    sort(H2Trees.aggregationnodes(testaggregateplan, level))
            end

            N = length(H2Trees.aggregationlevels(trialaggregatetranslateplan))
            for i in eachindex(testdisaggregationplan.nodes)
                tfd = H2Trees.translatingnodes(testdisaggregationplan)[i]
                tfa = H2Trees.receivingnodes(trialaggregatetranslateplan)[N - i + 1]

                for k in keys(tfd)
                    @test sort(tfd[k]) == sort(tfa[k])
                end
            end

            N = length(H2Trees.aggregationlevels(testaggregatetranslateplan))
            for i in eachindex(trialdisaggregationplan.nodes)
                tfd = H2Trees.translatingnodes(trialdisaggregationplan)[i]
                tfa = H2Trees.receivingnodes(testaggregatetranslateplan)[N - i + 1]

                for k in keys(tfd)
                    @test sort(tfd[k]) == sort(tfa[k])
                end
            end
        end
    end
end
