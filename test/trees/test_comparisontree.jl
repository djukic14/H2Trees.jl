using Test
using CompScienceMeshes
using H2Trees
@testset "ComparisonTwoNTree" begin
    meshes =
        [
            "cuboid",
            # "cuboid2",
            # "cuboid3",
            # "cuboid4",
            "multiplerects",
            # "multiplerects2",
            # "multiplerects3",
            "sphere",
            # "sphere2",
            # "sphere3",
            # "sphere4",
            # "sphere5",
            # "sphere6",
            # "sphere7",
            # "sphere8",
            "spherewithcenter",
            # "spherewithcenter2",
            # "spherewithcenter3",
            # "spherewithcenter4",
            # "spherewithcenter5",
            # "spherewithcenter6",
            # "spherewithcenter7",
            # "spherewithcenter8",
            # "spherewithcenter9",
            # "spherewithcenter10",
            # "spherewithcenter11",
            # "spherewithcenter12",
            # "spherewithcenter13",
            # "spherewithcenter14",
            # "spherewithcenter15",
            # "spherewithcenter16",
            "twospheres",
            # "twospheres2",
            # "twospheres3",
            # "twospheres4",
        ] .* ".in"

    ms = [
        CompScienceMeshes.readmesh(joinpath(pkgdir(H2Trees), "test", "assets", "in", m)) for
        m in meshes
    ]

    for minlevel in 1:3
        for root in 1:2
            for (i, m) in enumerate(ms)
                rootcenter, rootsize = H2Trees.boundingbox(vertices(m))

                tree = H2Trees.comparisonTwoNTree(vertices(m), root, rootsize, minlevel)

                @test H2Trees.checkbalancedtree(tree)
                @test length(H2Trees.leaves(tree)) == length(vertices(m))
                @test H2Trees.levels(tree)[1] == minlevel
                @test H2Trees.root(tree) == root

                valinnodes = [Int[] for i in 1:numvertices(m)]
                for leaf in H2Trees.leaves(tree)
                    @test length(H2Trees.values(tree, leaf)) == 1
                    @test H2Trees.isin(
                        tree, leaf, vertices(m)[H2Trees.values(tree, leaf)[1]]
                    )
                    for val in H2Trees.values(tree, leaf)
                        push!(valinnodes[val], leaf)
                    end
                end

                @test all(x -> length(x) == 1, valinnodes)

                for level in H2Trees.levels(tree)
                    for point in vertices(m)
                        node = H2Trees.locatepoint(tree, point, level)
                        @test H2Trees.isin(tree, node, point)
                        @test H2Trees.level(tree, node) == level
                    end
                end
                @test_throws ErrorException H2Trees.locatepoint(
                    tree, vertices(m)[1], H2Trees.levels(tree)[end] + 1
                )
            end
        end
    end
end
