using CompScienceMeshes, BEAST
using H2Trees
using Test
using SparseArrays
using ParallelKMeans
using LinearAlgebra
using PlotlyJS
@testset "KMeansTree" begin
    λ = 1.0

    ms = [
        CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere7.in")
        ),
        CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "cuboid3.in")
        ),
        CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "twospheres3.in")
        ),
        CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "multiplerects2.in")
        ),
    ]

    for i in eachindex(ms)
        tree = H2Trees.KMeansTree(vertices(ms[i]), 10; minvalues=100, n_threads=1)

        @test H2Trees.testwellseparatedness(tree)

        for leaf in H2Trees.leaves(tree)
            for point in H2Trees.values(tree, leaf)
                @test H2Trees.isin(tree, leaf, vertices(ms[i])[point])
            end
        end

        @test_nowarn println(tree)
        @test_nowarn display(tree)
        @test_nowarn show(tree)
    end

    for i in eachindex(ms)
        tree = H2Trees.KMeansTree(
            vertices(ms[i]),
            10;
            minvalues=100,
            n_threads=1,
            updateradii=H2Trees.noboundingsphereupdate,
        )

        for leaf in H2Trees.leaves(tree)
            for point in H2Trees.values(tree, leaf)
                @test H2Trees.isin(tree, leaf, vertices(ms[i])[point])
            end
        end

        @test_nowarn println(tree)
        @test_nowarn display(tree)
        @test_nowarn show(tree)
    end

    for i in eachindex(ms)
        tree = H2Trees.KMeansTree(
            vertices(ms[i]),
            10;
            minvalues=100,
            n_threads=1,
            updateradii=H2Trees.unsafemaxradiusboundingsphere,
        )

        for leaf in H2Trees.leaves(tree)
            for point in H2Trees.values(tree, leaf)
                @test H2Trees.isin(tree, leaf, vertices(ms[i])[point])
            end
        end

        for node in H2Trees.DepthFirstIterator(tree)
            radius = H2Trees.radius(tree, node)
            for child in H2Trees.children(tree, node)
                @test radius >= H2Trees.radius(tree, child)
            end
        end

        @test_nowarn println(tree)
        @test_nowarn display(tree)
        @test_nowarn show(tree)
    end
end
