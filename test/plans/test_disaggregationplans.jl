using Test
using BEAST, CompScienceMeshes
using StaticArrays
using LinearAlgebra
using H2Trees

@testset "Disaggregate Translate Plan" begin
    λ = 1.0
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere6.in")
    )
    X = raviartthomas(m)
    tree = buildtree(X; builder=TwoNTreeBuilder(; minhalfsize=λ / 10, minvalues=0))

    disaggregationplan = H2Trees.DisaggregateTranslatePlan(
        tree, H2Trees.TranslatingNodesIterator(; isnear=H2Trees.isnear())(tree)
    )
    # `plan[node, level]` yields the nodes translating into `node` at `level`, or an empty vector
    # when `node` receives nothing there.
    #
    # This used to be pinned as `disaggregationplan[10, 3] == Int[]`. That relied on node 10 not
    # being a level-3 node under the old depth-first numbering; node ids are now level-major and
    # Hilbert-ordered (level 3 starts exactly at node 10 here, and every level-3 node receives),
    # so a bare id encodes the layout rather than the behaviour under test. The two invariants it
    # was standing in for are asserted directly instead: layout-independent, and strictly
    # stronger since they cover every node rather than one.
    for level in H2Trees.levels(tree)
        receiving = Set(H2Trees.receivingnodes(disaggregationplan, level))
        for node in H2Trees.nodesatlevel(tree, level)
            if node in receiving
                @test !isempty(disaggregationplan[node, level])
            else
                @test isempty(disaggregationplan[node, level])
            end
            # A node asked about at a level that is not its own is never registered there.
            for otherlevel in H2Trees.levels(tree)
                otherlevel == level && continue
                @test isempty(disaggregationplan[node, otherlevel])
            end
        end
    end

    for node in H2Trees.DepthFirstIterator(tree, 1)
        nodehastobevisited = false

        if H2Trees.istranslatingnode(tree, node)
            nodehastobevisited = true
        end

        for parent in H2Trees.ParentUpwardsIterator(tree, node)
            if H2Trees.istranslatingnode(tree, parent)
                nodehastobevisited = true
                break
            end
        end

        if nodehastobevisited
            @test node in
                H2Trees.disaggregationnodes(disaggregationplan, H2Trees.level(tree, node))
            @test H2Trees.isdisaggregationnode(disaggregationplan, node)
        else
            @test node ∉ collect(
                Iterators.flatten(H2Trees.disaggregationnodes(disaggregationplan))
            )
            @test !H2Trees.isdisaggregationnode(disaggregationplan, node)
        end
    end

    for level in H2Trees.disaggregationlevels(disaggregationplan)
        for receivingnode in H2Trees.receivingnodes(disaggregationplan, level)
            @test H2Trees.translatingnodes(disaggregationplan, receivingnode, level) ==
                collect(H2Trees.TranslatingNodesIterator(tree, receivingnode))
        end
    end

    @test H2Trees.mintranslationlevel(disaggregationplan) == 3

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
            tree = buildtree(
                X,
                Y;
                builder=BlockTreeBuilder(;
                    test=TwoNTreeBuilder(; minhalfsize=λ / 10, minvalues=0),
                    trial=TwoNTreeBuilder(; minhalfsize=λ / 10, minvalues=0),
                ),
            )

            testtree = H2Trees.testtree(tree)
            trialtree = H2Trees.trialtree(tree)

            @test_throws ArgumentError H2Trees.DisaggregateTranslatePlan(
                tree, H2Trees.TranslatingNodesIterator
            )

            trialdisaggregationplan = H2Trees.DisaggregateTranslatePlan(
                trialtree, testtree, H2Trees.TranslatingNodesIterator
            )
            testdisaggregationplan = H2Trees.DisaggregateTranslatePlan(
                testtree, trialtree, H2Trees.TranslatingNodesIterator
            )

            for trialnode in H2Trees.DepthFirstIterator(trialtree, 1)
                nodehastobevisited = false

                if H2Trees.istranslatingnode(testtree, trialtree, trialnode)
                    nodehastobevisited = true
                end

                for parent in H2Trees.ParentUpwardsIterator(trialtree, trialnode)
                    if H2Trees.istranslatingnode(testtree, trialtree, parent)
                        nodehastobevisited = true
                        break
                    end
                end

                if nodehastobevisited
                    @test trialnode in H2Trees.disaggregationnodes(
                        trialdisaggregationplan, H2Trees.level(trialtree, trialnode)
                    )
                    @test H2Trees.isdisaggregationnode(trialdisaggregationplan, trialnode)
                else
                    @test trialnode ∉ collect(
                        Iterators.flatten(
                            H2Trees.disaggregationnodes(trialdisaggregationplan)
                        ),
                    )

                    @test !H2Trees.isdisaggregationnode(trialdisaggregationplan, trialnode)
                end
            end

            for level in H2Trees.disaggregationlevels(trialdisaggregationplan)
                for receivingnode in H2Trees.receivingnodes(trialdisaggregationplan, level)
                    @test trialdisaggregationplan[receivingnode, level] == collect(
                        H2Trees.TranslatingNodesIterator(testtree, trialtree, receivingnode)
                    )

                    @test trialdisaggregationplan[receivingnode, level] ==
                        H2Trees.translatingnodes(
                        trialdisaggregationplan, receivingnode, level
                    )

                    @test receivingnode in 1:H2Trees.numberofnodes(trialtree)
                end
            end

            for level in H2Trees.disaggregationlevels(testdisaggregationplan)
                for testnode in H2Trees.disaggregationnodes(testdisaggregationplan, level)
                    nodehastobevisited = false

                    if H2Trees.istranslatingnode(trialtree, testtree, testnode)
                        nodehastobevisited = true
                    end

                    for parent in H2Trees.ParentUpwardsIterator(testtree, testnode)
                        if H2Trees.istranslatingnode(trialtree, testtree, parent)
                            nodehastobevisited = true
                            break
                        end
                    end

                    if nodehastobevisited
                        @test testnode in H2Trees.disaggregationnodes(
                            testdisaggregationplan, H2Trees.level(testtree, testnode)
                        )
                    else
                        @test testnode ∉ collect(
                            Iterators.flatten(
                                H2Trees.disaggregationnodes(testdisaggregationplan)
                            ),
                        )
                    end

                    @test H2Trees.isdisaggregationnode(testdisaggregationplan, testnode)
                end

                for level in H2Trees.disaggregationlevels(testdisaggregationplan)
                    for receivingnode in
                        H2Trees.receivingnodes(testdisaggregationplan, level)
                        @test testdisaggregationplan[receivingnode, level] == collect(
                            H2Trees.TranslatingNodesIterator(
                                trialtree, testtree, receivingnode
                            ),
                        )

                        @test testdisaggregationplan[receivingnode, level] ==
                            H2Trees.translatingnodes(
                            testdisaggregationplan, receivingnode, level
                        )

                        @test receivingnode in 1:H2Trees.numberofnodes(testtree)
                    end
                end
            end
        end
    end
end

@testset "Disaggregate Plan" begin
    λ = 1.0
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere6.in")
    )
    X = raviartthomas(m)
    tree = buildtree(
        X;
        builder=TwoNTreeBuilder(;
            minhalfsize=λ / 10, minvalues=0, protrusion=NoProtrusionCheck()
        ),
    )

    TFIterator = H2Trees.TranslatingNodesIterator(; isnear=H2Trees.isnear())(tree)
    aggregatetranslateplan = H2Trees.AggregateTranslatePlan(tree, TFIterator)
    disaggregatetranslateplan = H2Trees.DisaggregateTranslatePlan(tree, TFIterator)
    compatibleplan = H2Trees.DisaggregateTranslatePlan(
        H2Trees.translatingnodes(disaggregatetranslateplan),
        H2Trees.nodes(disaggregatetranslateplan),
        H2Trees.levels(disaggregatetranslateplan),
        disaggregatetranslateplan.isdisaggregationnode,
        H2Trees.rootoffset(disaggregatetranslateplan),
        H2Trees.tree(disaggregatetranslateplan),
    )
    @test compatibleplan.receivingnodes_by_level ==
        disaggregatetranslateplan.receivingnodes_by_level
    refreshlevel = first(H2Trees.disaggregationlevels(compatibleplan))
    refreshlevelid = H2Trees.leveltolevelid(compatibleplan, refreshlevel)
    refreshdict = H2Trees.translatingnodes(compatibleplan)[refreshlevelid]
    refreshnode = maximum(collect(keys(refreshdict))) + 1
    refreshdict[refreshnode] = Int[]
    @test refreshnode ∉ H2Trees.receivingnodes(compatibleplan, refreshlevel)
    @test H2Trees.refreshreceivingnodes!(compatibleplan) === compatibleplan
    @test refreshnode in H2Trees.receivingnodes(compatibleplan, refreshlevel)
    delete!(refreshdict, refreshnode)
    H2Trees.refreshreceivingnodes!(compatibleplan)
    @test refreshnode ∉ H2Trees.receivingnodes(compatibleplan, refreshlevel)

    istranslatingfunctor = H2Trees.istranslatingnode(;
        TranslatingNodesIterator=H2Trees.TranslatingNodesIterator(; isnear=H2Trees.isnear())
    )(
        tree
    )

    aggregateplan = H2Trees.AggregatePlan(tree, istranslatingfunctor)
    disaggregateplan = H2Trees.DisaggregatePlan(tree, istranslatingfunctor)

    @test H2Trees.disaggregationlevels(disaggregateplan) ==
        H2Trees.disaggregationlevels(disaggregatetranslateplan)
    @test H2Trees.disaggregationlevels(disaggregateplan) ==
        reverse(H2Trees.aggregationlevels(aggregateplan))

    for level in H2Trees.disaggregationlevels(disaggregateplan)
        @test sort(H2Trees.disaggregationnodes(disaggregateplan, level)) ==
            sort(H2Trees.disaggregationnodes(disaggregatetranslateplan, level))
        @test sort(H2Trees.disaggregationnodes(disaggregateplan, level)) ==
            sort(H2Trees.aggregationnodes(aggregateplan, level))

        for node in H2Trees.disaggregationnodes(disaggregateplan, level)
            @test H2Trees.storenode(disaggregateplan, node) ==
                H2Trees.storenode(aggregateplan, node)
        end
    end

    @test aggregateplan.tree ==
        disaggregateplan.tree ==
        tree ==
        aggregatetranslateplan.tree ==
        disaggregatetranslateplan.tree

    @test H2Trees.rootoffset(disaggregateplan) ==
        H2Trees.rootoffset(disaggregatetranslateplan) ==
        H2Trees.rootoffset(aggregateplan) ==
        H2Trees.rootoffset(aggregatetranslateplan) ==
        H2Trees.root(tree) - 1

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
            tree = buildtree(
                X,
                Y;
                builder=BlockTreeBuilder(;
                    test=TwoNTreeBuilder(;
                        minhalfsize=λ / 10, minvalues=0, protrusion=NoProtrusionCheck()
                    ),
                    trial=TwoNTreeBuilder(;
                        minhalfsize=λ / 10, minvalues=0, protrusion=NoProtrusionCheck()
                    ),
                ),
            )
            TFIterator = H2Trees.TranslatingNodesIterator(; isnear=H2Trees.isnear())
            aggregatenode = H2Trees.istranslatingnode(; TranslatingNodesIterator=TFIterator)

            testtree = H2Trees.testtree(tree)
            trialtree = H2Trees.trialtree(tree)

            @test_throws ArgumentError H2Trees.AggregateTranslatePlan(
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

            trialdisaggregatetranslateplan = H2Trees.DisaggregateTranslatePlan(
                trialtree, testtree, H2Trees.TranslatingNodesIterator
            )
            testdisaggregatetranslateplan = H2Trees.DisaggregateTranslatePlan(
                testtree, trialtree, H2Trees.TranslatingNodesIterator
            )

            trialdisaggregateplan = H2Trees.DisaggregatePlan(
                trialtree,
                H2Trees.PetrovAggregationFunctor(aggregatenode, tree, testtree, trialtree),
            )

            testdisaggregateplan = H2Trees.DisaggregatePlan(
                testtree,
                H2Trees.PetrovAggregationFunctor(aggregatenode, tree, trialtree, testtree),
            )

            @test_throws ArgumentError H2Trees.DisaggregatePlan(
                tree, H2Trees.AggregateAllNodesFunctor()
            )

            @test H2Trees.disaggregationlevels(trialdisaggregateplan) ==
                H2Trees.disaggregationlevels(trialdisaggregatetranslateplan) ==
                reverse(H2Trees.aggregationlevels(testaggregateplan))

            @test H2Trees.disaggregationlevels(testdisaggregateplan) ==
                H2Trees.disaggregationlevels(testdisaggregatetranslateplan) ==
                reverse(H2Trees.aggregationlevels(trialaggregateplan))

            for level in H2Trees.disaggregationlevels(trialdisaggregateplan)
                @test sort(H2Trees.disaggregationnodes(trialdisaggregateplan, level)) ==
                    sort(
                    H2Trees.disaggregationnodes(trialdisaggregatetranslateplan, level)
                )
                @test sort(H2Trees.disaggregationnodes(trialdisaggregateplan, level)) ==
                    sort(H2Trees.aggregationnodes(trialaggregateplan, level))
            end

            for level in H2Trees.disaggregationlevels(testdisaggregateplan)
                @test sort(H2Trees.disaggregationnodes(testdisaggregateplan, level)) ==
                    sort(
                    H2Trees.disaggregationnodes(testdisaggregatetranslateplan, level)
                )
                @test sort(H2Trees.disaggregationnodes(testdisaggregateplan, level)) ==
                    sort(H2Trees.aggregationnodes(testaggregateplan, level))
            end

            @test trialaggregateplan.tree ==
                trialaggregatetranslateplan.tree ==
                trialdisaggregateplan.tree ==
                trialdisaggregatetranslateplan.tree ==
                trialtree

            @test testaggregateplan.tree ==
                testaggregatetranslateplan.tree ==
                testdisaggregateplan.tree ==
                testdisaggregatetranslateplan.tree ==
                testtree

            @test H2Trees.rootoffset(trialdisaggregateplan) ==
                H2Trees.rootoffset(trialdisaggregatetranslateplan) ==
                H2Trees.rootoffset(trialaggregateplan) ==
                H2Trees.rootoffset(trialaggregatetranslateplan) ==
                H2Trees.root(trialtree) - 1

            @test H2Trees.rootoffset(testdisaggregateplan) ==
                H2Trees.rootoffset(testdisaggregatetranslateplan) ==
                H2Trees.rootoffset(testaggregateplan) ==
                H2Trees.rootoffset(testaggregatetranslateplan) ==
                H2Trees.root(testtree) - 1
        end
    end
end
