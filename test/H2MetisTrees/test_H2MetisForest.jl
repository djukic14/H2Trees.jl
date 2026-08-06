
using Test
using BEAST, CompScienceMeshes
using Metis, Graphs
using H2Trees

@testset "Metis Forest" begin
    meshes =
        ["cuboid", "multiplerects", "sphere", "spherewithcenter", "twospheres"] .* ".in"

    ms = [
        CompScienceMeshes.readmesh(joinpath(pkgdir(H2Trees), "test", "assets", "in", m)) for
        m in meshes
    ]

    m = ms[1]
    _levelratios = []
    for (i, m) in enumerate(ms)
        # X = lagrangecxd0(m)
        for (j, X) in enumerate([lagrangecxd0(m), duallagrangecxd0(m), raviartthomas(m)])
            println("Testing mesh ($i / $(length(ms))) with function space ($j / 3)")
            g, w = H2Trees.adjacencygraph(X)

            forest = H2Trees.MetisForest(
                X;
                builder=MetisForestBuilder(;
                    treebuilder=MetisTreeBuilder(; numdivisions=4)
                ),
            )

            # Builder-only entry points over the (points, graph, weights) form.
            builtforest = H2Trees.buildforest(
                BEAST.positions(X),
                g,
                w;
                builder=MetisForestBuilder(;
                    treebuilder=MetisTreeBuilder(; numdivisions=4)
                ),
            )
            @test length(builtforest) == length(forest)
            spaceforest = H2Trees.buildforest(
                X;
                graphweights=(g, w),
                builder=MetisForestBuilder(;
                    treebuilder=MetisTreeBuilder(; numdivisions=4)
                ),
            )
            @test length(spaceforest) == length(forest)
            builttree = H2Trees.buildtree(
                BEAST.positions(X), g, w; builder=MetisTreeBuilder(; numdivisions=4)
            )
            @test builttree isa H2Trees.BoundingBallTree
            spacetree = H2Trees.buildtree(
                X; graphweights=(g, w), builder=MetisTreeBuilder(; numdivisions=4)
            )
            directspacetree = H2Trees.MetisTree(
                X; graphweights=(g, w), builder=MetisTreeBuilder(; numdivisions=4)
            )
            @test spacetree isa H2Trees.BoundingBallTree
            @test H2Trees.numberofnodes(spacetree) == H2Trees.numberofnodes(builttree)
            @test H2Trees.numberofnodes(directspacetree) == H2Trees.numberofnodes(spacetree)

            vals = Int[]

            @test length(forest) == length(forest.trees)
            eltype(forest) == typeof(forest.trees[1])
            @test forest[1] == forest.trees[1]

            for tree in forest
                for node in H2Trees.DepthFirstIterator(tree)
                    @test length(collect(H2Trees.children(tree, node))) <= 4

                    subgraph, _ = induced_subgraph(g, H2Trees.values(tree, node))
                    @test is_connected(subgraph)
                end

                for leaf in H2Trees.leaves(tree)
                    @test length(H2Trees.values(tree, leaf)) <= 4
                    append!(vals, H2Trees.values(tree, leaf))
                end
            end
            sort!(vals)
            @test vals == 1:numfunctions(X)

            compact = sprint(show, forest)
            @test occursin("Forest(trees=", compact)
            @test occursin("nodes=", compact)
            @test occursin("leaves=", compact)
            @test occursin("levels=", compact)
            @test !occursin('\n', compact)

            detailed = sprint(show, MIME"text/plain"(), forest)
            @test occursin("Forest with $(length(forest)) tree(s)", detailed)
            @test occursin("[1] ", detailed)
        end
    end
end

@testset "H2METISTrees internals" begin
    H2MetisTrees = Base.get_extension(H2Trees, :H2MetisTrees)

    meshes =
        ["cuboid", "multiplerects", "sphere", "spherewithcenter", "twospheres"] .* ".in"

    ms = [
        CompScienceMeshes.readmesh(joinpath(pkgdir(H2Trees), "test", "assets", "in", m)) for
        m in meshes
    ]

    for m in ms
        X = lagrangecxd0(m)

        g, w = H2Trees.adjacencygraph(X)

        res = H2MetisTrees.metispartition(g, w, 4; alg=:RECURSIVE)
        @test length(res) == numfunctions(X)
        @test maximum(res) <= 4

        @test H2MetisTrees.computemetisweights(Int, zeros(numfunctions(X))) ==
            ones(Int, numfunctions(X))
    end
end
