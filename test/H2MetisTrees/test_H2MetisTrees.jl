using Test
using BEAST, CompScienceMeshes
using Metis, Graphs
using H2Trees
using Test
using BEAST, CompScienceMeshes
using Metis, Graphs
using H2Trees
@testset "Metis trees" begin
    meshes =
        ["cuboid", "multiplerects", "sphere", "spherewithcenter", "twospheres"] .* ".in"

    ms = [
        CompScienceMeshes.readmesh(joinpath(pkgdir(H2Trees), "test", "assets", "in", m)) for
        m in meshes
    ]

    m = ms[1]
    _levelratios = []
    for (i, m) in enumerate(ms)
        for (j, X) in enumerate([lagrangecxd0(m), duallagrangecxd0(m), raviartthomas(m)])
            for split in [true, false]
                println("Testing mesh ($i / $(length(ms))) with function space ($j / 3)")
                tree = H2Trees.MetisTree(
                    X;
                    builder=MetisTreeBuilder(;
                        numdivisions=4, splitunconnectedpartitions=split
                    ),
                )

                for node in H2Trees.DepthFirstIterator(tree)
                    if !split
                        @test length(collect(H2Trees.children(tree, node))) <= 4
                    end
                end

                vals = Int[]
                for leaf in H2Trees.leaves(tree)
                    @test length(H2Trees.values(tree, leaf)) <= 4
                    append!(vals, H2Trees.values(tree, leaf))
                end

                sort!(vals)
                @test vals == 1:numfunctions(X)

                elements, ad, _ = assemblydata(X; onlyactives=false)
                functionsperelement = [length(ad[i]) for i in 1:length(elements)]

                area = sum(
                    CompScienceMeshes.volume(chart(X.geo, i)) for i in 1:numcells(X.geo)
                )
                areas = [
                    CompScienceMeshes.volume(chart(X.geo, i)) / functionsperelement[i] for
                    i in 1:numcells(X.geo)
                ]

                nodearea = zeros(H2Trees.numberofnodes(tree))
                for node in 1:H2Trees.numberofnodes(tree)
                    for val in H2Trees.values(tree, node)
                        for fn in X.fns[val]
                            nodearea[node] += areas[fn.cellid]
                        end
                    end
                end
                levelareas = [
                    (nodearea[H2Trees.LevelIterator(tree, level)]) for
                    level in H2Trees.levels(tree)
                ]

                leaflevel = minimum(
                    sort!(unique(H2Trees.level.(Ref(tree), H2Trees.leaves(tree))))
                )

                maxminlevelareas = [
                    (maximum(levelareas[i]), minimum(levelareas[i])) for
                    i in eachindex(levelareas)
                ]
                levelratios = [
                    maxminlevelareas[i][1] / maxminlevelareas[i][2] for
                    i in eachindex(maxminlevelareas)
                ]
                println("Max min level area ratios: ", levelratios)
                push!(_levelratios, levelratios)

                @test maximum(levelratios[1:(leaflevel - 1)]) < 2
                b = boundary(m)

                if isempty(b)
                    for i in eachindex(levelareas)
                        i <= leaflevel && @test sum(levelareas[i]) ≈ area
                    end
                end

                if split
                    g, w = H2Trees.adjacencygraph(X)
                    for node in H2Trees.DepthFirstIterator(tree)
                        node == H2Trees.root(tree) && continue
                        subgraph, _ = induced_subgraph(g, H2Trees.values(tree, node))
                        @test is_connected(subgraph)
                    end
                end
            end
        end
    end

    @testset "balanceleaves! on realistic Metis geometries" begin
        sawunbalanced = false
        for m in ms
            X = raviartthomas(m)
            tree = H2Trees.MetisTree(
                X; builder=MetisTreeBuilder(; numdivisions=4, minvalues=8)
            )

            beforelevels = sort(unique(H2Trees.level.(Ref(tree), H2Trees.leaves(tree))))
            sawunbalanced |= length(beforelevels) > 1
            beforevalues = sort!(
                reduce(vcat, H2Trees.values.(Ref(tree), H2Trees.leaves(tree)))
            )

            @test H2Trees.balanceleaves!(tree) === tree
            @test H2Trees.checkbalancedtree(tree)
            @test length(unique(H2Trees.level.(Ref(tree), H2Trees.leaves(tree)))) == 1
            @test sort!(reduce(vcat, H2Trees.values.(Ref(tree), H2Trees.leaves(tree)))) ==
                beforevalues
            @test H2Trees.nodesatlevel(tree) == H2Trees.treeindex(tree).nodes_by_level
            @test H2Trees.treeindex(tree).leaves == H2Trees.leaves(tree)
            @test H2Trees.depthfirstnodes(tree) ==
                collect(Int, H2Trees.DepthFirstIterator(tree))

            for node in H2Trees.DepthFirstIterator(tree)
                if !H2Trees.isleaf(tree, node)
                    @test isempty(H2Trees.values(H2Trees.data(tree, node)))
                end
            end
            for leaf in H2Trees.leaves(tree)
                for value in H2Trees.values(tree, leaf)
                    @test H2Trees.isin(tree, leaf, BEAST.positions(X)[value])
                end
            end
        end
        @test sawunbalanced
    end
end
