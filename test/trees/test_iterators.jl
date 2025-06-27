using Test
using CompScienceMeshes
using StaticArrays
using H2Trees
using BEAST

@testset "Depth First Iterator" begin
    λ = 1.0

    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "twospheres.in")
    )
    X = raviartthomas(m)

    minlevels = 1:10
    roots = 1:3
    minvalues = [0, 11, 20]

    for root in roots
        for minlevel in minlevels
            for minvalues in minvalues
                tree = TwoNTree(X, λ / 5; minlevel=minlevel, root=root, minvalues=minvalues)

                valuesatnodes = H2Trees.valuesatnodes(tree)
                @test length(valuesatnodes) == numfunctions(X)
                for (functionid, value) in enumerate(valuesatnodes)
                    @test length(value) == 1
                    @test functionid in H2Trees.values(tree, value[1])
                end
                nodesatvalues = H2Trees.nodesatvalues(tree)
                for (key, value) in nodesatvalues
                    @test length(key) == 1
                    key = key[1]
                    @test sort(value) == sort(H2Trees.values(tree, key))
                end

                @test H2Trees.minimumlevel(tree) == minlevel
                @test H2Trees.levels(tree)[1] == minlevel

                @test sort!(collect(H2Trees.DepthFirstIterator(tree, root))) ==
                    Array(root:(length(tree.nodes) + root - 1))
                @test typeof(collect(H2Trees.DepthFirstIterator(tree, root))) == Vector{Int}

                for i in H2Trees.LevelIterator(tree, minlevel + 1)
                    @test H2Trees.level(tree, i) == minlevel + 1
                    collect(H2Trees.ParentUpwardsIterator(tree, i)) == [1]
                    @test typeof(collect(H2Trees.ParentUpwardsIterator(tree, i))) ==
                        Vector{Int}
                end
            end
        end
    end
end

@testset "Parent Upwards" begin
    λ = 1.0
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere.in")
    )
    X = raviartthomas(m)

    root = 3
    for minvalues in [0, 10]
        tree = TwoNTree(X, λ / 10; minlevel=root, minvalues=minvalues)

        valuesatnodes = H2Trees.valuesatnodes(tree)
        @test length(valuesatnodes) == numfunctions(X)
        for (functionid, value) in enumerate(valuesatnodes)
            @test length(value) == 1
            @test functionid in H2Trees.values(tree, value[1])
        end
        nodesatvalues = H2Trees.nodesatvalues(tree)
        for (key, value) in nodesatvalues
            @test length(key) == 1
            key = key[1]
            @test sort(value) == sort(H2Trees.values(tree, key))
        end

        for level in H2Trees.levels(tree)
            for node in H2Trees.LevelIterator(tree, level)
                nodes = collect(H2Trees.ParentUpwardsIterator(tree, node))

                level == 2 && @test nodes == []

                for i in eachindex(nodes)
                    if i < length(nodes)
                        @test nodes[i + 1] == tree(nodes[i]).parent
                    end
                end
            end
        end
    end
end

@testset "Near/FarNodes TwoNTree" begin
    λ = 1.0
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere.in")
    )

    X = raviartthomas(m)

    root = 3
    for minvalues in [0, 10]
        tree = TwoNTree(X, λ / 10; minlevel=2, root=root, minvalues=minvalues)

        valuesatnodes = H2Trees.valuesatnodes(tree)
        @test length(valuesatnodes) == numfunctions(X)
        for (functionid, value) in enumerate(valuesatnodes)
            @test length(value) == 1
            @test functionid in H2Trees.values(tree, value[1])
        end
        nodesatvalues = H2Trees.nodesatvalues(tree)
        for (key, value) in nodesatvalues
            @test length(key) == 1
            key = key[1]
            @test sort(value) == sort(H2Trees.values(tree, key))
        end

        nearfunctor = H2Trees.isnear()(tree)
        farfunctor = H2Trees.isfar(nearfunctor)
        for node in H2Trees.DepthFirstIterator(tree)
            for nearnode in H2Trees.NearNodeIterator(tree, node)
                @test nearfunctor(tree, node, nearnode)
                @test !farfunctor(tree, node, nearnode)
            end
        end

        for centernode in eachindex(tree.nodes)
            centernode = centernode + root - 1
            for node in H2Trees.TranslatingNodesIterator(tree, centernode)
                @test H2Trees.level(tree, node) == H2Trees.level(tree, centernode)
                @test !H2Trees.isnear(tree, node, centernode)

                parentsnodes = collect(H2Trees.ParentUpwardsIterator(tree, node))
                parentscenternode = collect(H2Trees.ParentUpwardsIterator(tree, centernode))

                for i in eachindex(parentsnodes)
                    @test H2Trees.level(tree, parentsnodes[i]) ==
                        H2Trees.level(tree, parentscenternode[i])
                    @test H2Trees.isnear(tree, parentsnodes[i], parentscenternode[i])
                end
            end

            translatingnodes = collect(H2Trees.TranslatingNodesIterator(tree, centernode))
            nottranslatingnodes = collect(
                H2Trees.NotTranslatingNodesIterator(tree, centernode)
            )

            @test sort([translatingnodes; nottranslatingnodes]) ==
                Array(H2Trees.LevelIterator(tree, H2Trees.level(tree, centernode)))

            for node in H2Trees.FarNodeIterator(tree, centernode)
                @test H2Trees.level(tree, node) == H2Trees.level(tree, centernode)
                @test !H2Trees.isnear(
                    tree, node, centernode, H2Trees.isTwoNTree(); additionalbufferboxes=0
                )
            end

            for node in H2Trees.NearNodeIterator(tree, centernode)
                @test H2Trees.level(tree, node) == H2Trees.level(tree, centernode)
                @test H2Trees.isnear(
                    tree, node, centernode, H2Trees.isTwoNTree(); additionalbufferboxes=0
                )
            end

            nearnodevalues = H2Trees.nearnodevalues(tree, centernode)
            farnodevalues = H2Trees.farnodevalues(tree, centernode)

            @test sort([nearnodevalues; farnodevalues]) == Array(1:numfunctions(X))

            @test typeof(collect(H2Trees.FarNodeIterator(tree, centernode))) == Vector{Int}

            @test typeof(collect(H2Trees.NearNodeIterator(tree, centernode))) == Vector{Int}
        end

        for additionalbufferboxes in [0, 1, 2]
            isfarfilter =
                (tree, node, testnode) ->
                    !H2Trees.isnear(
                        tree,
                        node,
                        testnode,
                        H2Trees.isTwoNTree();
                        additionalbufferboxes=additionalbufferboxes,
                    )
            isnearfilter =
                (tree, node, testnode) -> H2Trees.isnear(
                    tree,
                    node,
                    testnode,
                    H2Trees.isTwoNTree();
                    additionalbufferboxes=additionalbufferboxes,
                )

            iswellseparatedfilter =
                (tree, node, testnode) -> H2Trees.iswellseparated(
                    tree, node, testnode, H2Trees.isTwoNTree(); isnear=isnearfilter
                )

            isnotwellseparatedfilter =
                (tree, node, testnode) ->
                    !H2Trees.iswellseparated(
                        tree, node, testnode, H2Trees.isTwoNTree(); isnear=isnearfilter
                    )

            for centernode in eachindex(tree.nodes)
                centernode = centernode + root - 1
                translatingnodes = collect(
                    H2Trees.TranslatingNodesIterator(
                        tree, centernode; iswellseparated=iswellseparatedfilter
                    ),
                )

                nottranslatingnodes = collect(
                    H2Trees.NotTranslatingNodesIterator(
                        tree, centernode; isnotwellseparated=isnotwellseparatedfilter
                    ),
                )

                @test sort([translatingnodes; nottranslatingnodes]) ==
                    Array(H2Trees.LevelIterator(tree, H2Trees.level(tree, centernode)))

                for testnode in H2Trees.TranslatingNodesIterator(
                    tree, centernode; iswellseparated=iswellseparatedfilter
                )
                    @test H2Trees.level(tree, testnode) == H2Trees.level(tree, centernode)

                    @test !H2Trees.isnear(
                        tree,
                        testnode,
                        centernode,
                        H2Trees.isTwoNTree();
                        additionalbufferboxes=additionalbufferboxes,
                    )

                    parentstestnode = collect(H2Trees.ParentUpwardsIterator(tree, testnode))
                    parentscenternode = collect(
                        H2Trees.ParentUpwardsIterator(tree, centernode)
                    )

                    for i in eachindex(parentstestnode)
                        @test H2Trees.level(tree, parentstestnode[i]) ==
                            H2Trees.level(tree, parentscenternode[i])
                        @test H2Trees.isnear(
                            tree,
                            parentstestnode[i],
                            parentscenternode[i],
                            H2Trees.isTwoNTree();
                            additionalbufferboxes=additionalbufferboxes,
                        )
                    end
                end
                for node in H2Trees.FarNodeIterator(tree, centernode; isfar=isfarfilter)
                    @test H2Trees.level(tree, node) == H2Trees.level(tree, centernode)
                    @test !H2Trees.isnear(
                        tree,
                        node,
                        centernode,
                        H2Trees.isTwoNTree();
                        additionalbufferboxes=additionalbufferboxes,
                    )
                end

                for node in H2Trees.NearNodeIterator(tree, centernode; isnear=isnearfilter)
                    @test H2Trees.level(tree, node) == H2Trees.level(tree, centernode)
                    @test H2Trees.isnear(
                        tree,
                        node,
                        centernode,
                        H2Trees.isTwoNTree();
                        additionalbufferboxes=additionalbufferboxes,
                    )
                end

                nearnodevalues = H2Trees.nearnodevalues(
                    tree, centernode; isnear=isnearfilter
                )
                farnodevalues = H2Trees.farnodevalues(tree, centernode; isfar=isfarfilter)

                @test sort([nearnodevalues; farnodevalues]) == Array(1:numfunctions(X))

                @test typeof(
                    collect(H2Trees.FarNodeIterator(tree, centernode; isfar=isfarfilter))
                ) == Vector{Int}

                @test typeof(
                    collect(H2Trees.NearNodeIterator(tree, centernode; isnear=isnearfilter))
                ) == Vector{Int}
            end
        end
    end
end

@testset "Near/FarNodes BlockTree" begin
    λ = 1.0
    mx = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter.in")
    )
    my = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter2.in")
    )

    X = raviartthomas(mx)
    Y = raviartthomas(my)

    tree = TwoNTree(X, Y, λ / 10)
    testtree = H2Trees.testtree(tree)
    trialtree = H2Trees.trialtree(tree)

    isnearfunctor = H2Trees.isnear()(tree)
    isfarfunctor = H2Trees.isfar(isnearfunctor)

    for trialnode in H2Trees.DepthFirstIterator(trialtree)
        for testnode in H2Trees.NearNodeIterator(testtree, trialtree, trialnode)
            @test isnearfunctor(testtree, trialtree, testnode, trialnode)
            @test !isfarfunctor(testtree, trialtree, testnode, trialnode)
        end
    end

    for tree in [testtree, trialtree]
        valuesatnodes = H2Trees.valuesatnodes(tree)
        for (functionid, value) in enumerate(valuesatnodes)
            @test length(value) == 1
            @test functionid in H2Trees.values(tree, value[1])
        end
        nodesatvalues = H2Trees.nodesatvalues(tree)
        for (key, value) in nodesatvalues
            @test length(key) == 1
            key = key[1]
            @test sort(value) == sort(H2Trees.values(tree, key))
        end
    end

    @test H2Trees.levels(testtree) == 2:5
    @test H2Trees.levels(trialtree) == 1:5

    @test_nowarn println(tree)
    @test_nowarn show(tree)
    @test_nowarn display(tree)

    for trialnode in eachindex(H2Trees.trialtree(tree).nodes)
        for testnode in H2Trees.TranslatingNodesIterator(tree, trialnode)
            @test H2Trees.level(testtree, testnode) == H2Trees.level(trialtree, trialnode)
            @test !H2Trees.isnear(
                testtree,
                trialtree,
                testnode,
                trialnode,
                H2Trees.isTwoNTree(),
                H2Trees.isTwoNTree();
            )
            parentstestnode = collect(H2Trees.ParentUpwardsIterator(testtree, testnode))
            parentstrialnode = collect(H2Trees.ParentUpwardsIterator(trialtree, trialnode))

            parentstestnode == [] && continue
            parentstrialnode == [] && continue
            for i in eachindex(parentstestnode)
                @test H2Trees.level(testtree, parentstestnode[i]) ==
                    H2Trees.level(trialtree, parentstrialnode[i])
                @test H2Trees.isnear(
                    testtree,
                    trialtree,
                    parentstestnode[i],
                    parentstrialnode[i],
                    H2Trees.isTwoNTree(),
                    H2Trees.isTwoNTree();
                )
            end
        end

        translatingnodes = collect(H2Trees.TranslatingNodesIterator(tree, trialnode))
        nottranslatingnodes = collect(H2Trees.NotTranslatingNodesIterator(tree, trialnode))

        @test sort([translatingnodes; nottranslatingnodes]) ==
            Array(H2Trees.LevelIterator(testtree, H2Trees.level(trialtree, trialnode)))

        @test collect(H2Trees.NotTranslatingNodesIterator(tree, trialnode)) ==
            collect(H2Trees.NotTranslatingNodesIterator(testtree, trialtree, trialnode))

        @test collect(H2Trees.NearNodeIterator(tree, trialnode)) ==
            collect(H2Trees.NearNodeIterator(testtree, trialtree, trialnode))
        for testnode in H2Trees.FarNodeIterator(tree, trialnode)
            @test H2Trees.level(testtree, testnode) == H2Trees.level(trialtree, trialnode)
            @test !H2Trees.isnear(
                testtree,
                trialtree,
                testnode,
                trialnode,
                H2Trees.isTwoNTree(),
                H2Trees.isTwoNTree();
            )
        end

        for testnode in H2Trees.NearNodeIterator(tree, trialnode)
            @test H2Trees.level(testtree, testnode) == H2Trees.level(trialtree, trialnode)
            @test H2Trees.isnear(
                testtree,
                trialtree,
                testnode,
                trialnode,
                H2Trees.isTwoNTree(),
                H2Trees.isTwoNTree();
            )
        end

        @test collect(H2Trees.FarNodeIterator(tree, trialnode)) ==
            collect(H2Trees.FarNodeIterator(testtree, trialtree, trialnode))

        @test collect(H2Trees.NearNodeIterator(tree, trialnode)) ==
            collect(H2Trees.NearNodeIterator(testtree, trialtree, trialnode))
    end

    for testnode in eachindex(H2Trees.testtree(tree).nodes)
        for trialnode in H2Trees.FarNodeIterator(trialtree, testtree, testnode)
            @test H2Trees.level(testtree, testnode) == H2Trees.level(trialtree, trialnode)
            @test !H2Trees.isnear(
                testtree,
                trialtree,
                testnode,
                trialnode,
                H2Trees.isTwoNTree(),
                H2Trees.isTwoNTree();
            )
        end

        for trialnode in H2Trees.NearNodeIterator(trialtree, testtree, testnode)
            @test H2Trees.level(testtree, testnode) == H2Trees.level(trialtree, trialnode)
            @test H2Trees.isnear(
                testtree,
                trialtree,
                testnode,
                trialnode,
                H2Trees.isTwoNTree(),
                H2Trees.isTwoNTree();
            )
        end
    end

    for additionalbufferboxes in [0, 1, 2]
        isfarfilter =
            (testtree, trialtree, trialnode, testnode) ->
                !H2Trees.isnear(
                    testtree,
                    trialtree,
                    trialnode,
                    testnode,
                    H2Trees.isTwoNTree(),
                    H2Trees.isTwoNTree();
                    additionalbufferboxes=additionalbufferboxes,
                )

        isnearfilter =
            (testtree, trialtree, trialnode, testnode) -> H2Trees.isnear(
                testtree,
                trialtree,
                trialnode,
                testnode,
                H2Trees.isTwoNTree(),
                H2Trees.isTwoNTree();
                additionalbufferboxes=additionalbufferboxes,
            )

        iswellseparatedfilter =
            (testtree, trialtree, testnode, trialnode) -> H2Trees.iswellseparated(
                testtree,
                trialtree,
                testnode,
                trialnode,
                H2Trees.isTwoNTree(),
                H2Trees.isTwoNTree();
                isnear=isnearfilter,
            )

        isnotwellseparatedfilter =
            (testtree, trialtree, testnode, trialnode) ->
                !H2Trees.iswellseparated(
                    testtree,
                    trialtree,
                    testnode,
                    trialnode,
                    H2Trees.isTwoNTree(),
                    H2Trees.isTwoNTree();
                    isnear=isnearfilter,
                )

        for trialnode in eachindex(H2Trees.trialtree(tree).nodes)
            translatingnodes = collect(
                H2Trees.TranslatingNodesIterator(
                    tree, trialnode; iswellseparated=iswellseparatedfilter
                ),
            )
            nottranslatingnodes = collect(
                H2Trees.NotTranslatingNodesIterator(
                    tree, trialnode; isnotwellseparated=isnotwellseparatedfilter
                ),
            )

            @test sort([translatingnodes; nottranslatingnodes]) == Array(
                H2Trees.LevelIterator(testtree, H2Trees.level(trialtree, trialnode))
            )

            @test collect(
                H2Trees.NotTranslatingNodesIterator(
                    tree, trialnode; isnotwellseparated=isnotwellseparatedfilter
                ),
            ) == collect(
                H2Trees.NotTranslatingNodesIterator(
                    testtree,
                    trialtree,
                    trialnode;
                    isnotwellseparated=isnotwellseparatedfilter,
                ),
            )

            for testnode in H2Trees.TranslatingNodesIterator(
                tree, trialnode; iswellseparated=iswellseparatedfilter
            )
                @test H2Trees.level(testtree, testnode) ==
                    H2Trees.level(trialtree, trialnode)
                @test !H2Trees.isnear(
                    testtree,
                    trialtree,
                    testnode,
                    trialnode,
                    H2Trees.isTwoNTree(),
                    H2Trees.isTwoNTree();
                    additionalbufferboxes=additionalbufferboxes,
                )

                parentstestnode = collect(H2Trees.ParentUpwardsIterator(testtree, testnode))
                parentstrialnode = collect(
                    H2Trees.ParentUpwardsIterator(trialtree, trialnode)
                )

                parentstestnode == [] && continue
                parentstrialnode == [] && continue
                for i in eachindex(parentstestnode)
                    @test H2Trees.level(testtree, parentstestnode[i]) ==
                        H2Trees.level(trialtree, parentstrialnode[i])
                    @test H2Trees.isnear(
                        testtree,
                        trialtree,
                        parentstestnode[i],
                        parentstrialnode[i],
                        H2Trees.isTwoNTree(),
                        H2Trees.isTwoNTree();
                        additionalbufferboxes=additionalbufferboxes,
                    )
                end
            end

            for testnode in H2Trees.FarNodeIterator(tree, trialnode; isfar=isfarfilter)
                @test H2Trees.level(testtree, testnode) ==
                    H2Trees.level(trialtree, trialnode)
                @test !H2Trees.isnear(
                    testtree,
                    trialtree,
                    testnode,
                    trialnode,
                    H2Trees.isTwoNTree(),
                    H2Trees.isTwoNTree();
                    additionalbufferboxes=additionalbufferboxes,
                )
            end

            for testnode in H2Trees.NearNodeIterator(tree, trialnode; isnear=isnearfilter)
                @test H2Trees.level(testtree, testnode) ==
                    H2Trees.level(trialtree, trialnode)
                @test H2Trees.isnear(
                    testtree,
                    trialtree,
                    testnode,
                    trialnode,
                    H2Trees.isTwoNTree(),
                    H2Trees.isTwoNTree();
                    additionalbufferboxes=additionalbufferboxes,
                )
            end

            @test typeof(
                collect(H2Trees.FarNodeIterator(tree, trialnode; isfar=isfarfilter))
            ) == Vector{Int}

            @test typeof(
                collect(H2Trees.NearNodeIterator(tree, trialnode; isnear=isnearfilter))
            ) == Vector{Int}

            @test collect(H2Trees.FarNodeIterator(tree, trialnode; isfar=isfarfilter)) ==
                collect(
                H2Trees.FarNodeIterator(testtree, trialtree, trialnode; isfar=isfarfilter)
            )

            @test collect(H2Trees.NearNodeIterator(tree, trialnode; isnear=isnearfilter)) ==
                collect(
                H2Trees.NearNodeIterator(
                    testtree, trialtree, trialnode; isnear=isnearfilter
                ),
            )

            H2Trees.level(trialtree, trialnode) ∉ H2Trees.levels(testtree) && continue
            nearnodevalues = H2Trees.nearnodevalues(
                testtree, trialtree, trialnode; isnear=isnearfilter
            )
            farnodevalues = H2Trees.farnodevalues(
                testtree, trialtree, trialnode; isfar=isfarfilter
            )

            @test sort([nearnodevalues; farnodevalues]) == Array(1:numfunctions(X))
        end

        for testnode in eachindex(H2Trees.testtree(tree).nodes)
            for trialnode in
                H2Trees.FarNodeIterator(trialtree, testtree, testnode; isfar=isfarfilter)
                @test H2Trees.level(testtree, testnode) ==
                    H2Trees.level(trialtree, trialnode)
                @test !H2Trees.isnear(
                    testtree,
                    trialtree,
                    testnode,
                    trialnode,
                    H2Trees.isTwoNTree(),
                    H2Trees.isTwoNTree();
                    additionalbufferboxes=additionalbufferboxes,
                )
            end

            for trialnode in
                H2Trees.NearNodeIterator(trialtree, testtree, testnode; isnear=isnearfilter)
                @test H2Trees.level(testtree, testnode) ==
                    H2Trees.level(trialtree, trialnode)
                @test H2Trees.isnear(
                    testtree,
                    trialtree,
                    testnode,
                    trialnode,
                    H2Trees.isTwoNTree(),
                    H2Trees.isTwoNTree();
                    additionalbufferboxes=additionalbufferboxes,
                )
            end

            @test typeof(
                collect(
                    H2Trees.FarNodeIterator(
                        trialtree, testtree, testnode; isfar=isfarfilter
                    ),
                ),
            ) == Vector{Int}

            @test typeof(
                collect(
                    H2Trees.NearNodeIterator(
                        trialtree, testtree, testnode; isnear=isnearfilter
                    ),
                ),
            ) == Vector{Int}

            H2Trees.level(testtree, testnode) ∉ H2Trees.levels(trialtree) && continue

            nearnodevalues = H2Trees.nearnodevalues(
                trialtree, testtree, testnode; isnear=isnearfilter
            )
            farnodevalues = H2Trees.farnodevalues(
                trialtree, testtree, testnode; isfar=isfarfilter
            )

            @test sort([nearnodevalues; farnodevalues]) == Array(1:numfunctions(Y))
        end
    end
end

@testset "Connections of basis functions" begin
    λ = 1.0
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere2.in")
    )
    X = raviartthomas(m)

    root = 4
    for minvalues in [0, 10]
        tree = TwoNTree(X, λ / 10; minlevel=2, root=root, minvalues=minvalues)

        valuesatnodes = H2Trees.valuesatnodes(tree)
        @test length(valuesatnodes) == numfunctions(X)
        for (functionid, value) in enumerate(valuesatnodes)
            @test length(value) == 1
            @test functionid in H2Trees.values(tree, value[1])
        end
        nodesatvalues = H2Trees.nodesatvalues(tree)
        for (key, value) in nodesatvalues
            @test length(key) == 1
            key = key[1]
            @test sort(value) == sort(H2Trees.values(tree, key))
        end

        @test H2Trees.testwellseparatedness(tree)

        @test H2Trees.WellSeparatedIterator(; iswellseparated=1).iswellseparated == 1
    end
end

@testset "Connections basis test functions" begin
    λ = 1.0
    sphereλ = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter3.in")
    )
    sphere2λ = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter4.in")
    )
    sphere2λdisplaced = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter5.in")
    )
    sphere2λdisplacedfar = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter6.in")
    )

    ms = [sphereλ, sphere2λ, sphere2λdisplaced, sphere2λdisplacedfar]

    for minvaluestest in [0, 5]
        for minvaluestrial in [0, 5]
            for mx in ms
                X = raviartthomas(mx)
                for my in ms
                    Y = raviartthomas(my)
                    tree = TwoNTree(
                        X,
                        Y,
                        λ / 20;
                        minvaluestest=minvaluestrial,
                        minvaluestrial=minvaluestrial,
                    )

                    testlevels = unique(
                        H2Trees.level.(
                            Ref(H2Trees.testtree(tree)),
                            H2Trees.leaves(H2Trees.testtree(tree)),
                        ),
                    )
                    triallevels = unique(
                        H2Trees.level.(
                            Ref(H2Trees.trialtree(tree)),
                            H2Trees.leaves(H2Trees.trialtree(tree)),
                        ),
                    )

                    for tree in [H2Trees.testtree(tree), H2Trees.trialtree(tree)]
                        valuesatnodes = H2Trees.valuesatnodes(tree)
                        for (functionid, value) in enumerate(valuesatnodes)
                            @test length(value) == 1
                            @test functionid in H2Trees.values(tree, value[1])
                        end
                        nodesatvalues = H2Trees.nodesatvalues(tree)
                        for (key, value) in nodesatvalues
                            @test length(key) == 1
                            key = key[1]
                            @test sort(value) == sort(H2Trees.values(tree, key))
                        end
                    end
                    @test H2Trees.testwellseparatedness(tree)
                end
            end
        end
    end
end
