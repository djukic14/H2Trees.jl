module TestH2ParallelKMeansTrees

using CompScienceMeshes, BEAST
using H2Trees
using Test
using SparseArrays
using ParallelKMeans
using LinearAlgebra
using PlotlyJS
using Random
using StaticArrays

include(joinpath(pkgdir(H2Trees), "test", "testutils.jl"))

@testset verbose = true "KMeansTree" begin
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
        tree = H2Trees.buildtree(
            vertices(ms[i]);
            builder=H2Trees.KMeansTreeBuilder(;
                numberofclusters=10,
                minvalues=100,
                splitterkwargs=(; n_threads=1, rng=Random.MersenneTwister(1234)),
            ),
        )

        @test TestingUtils.testwellseparatedness(tree)

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
        tree = H2Trees.buildtree(
            vertices(ms[i]);
            builder=H2Trees.KMeansTreeBuilder(;
                numberofclusters=10,
                minvalues=100,
                splitterkwargs=(; n_threads=1, rng=Random.MersenneTwister(1234)),
                updateradii=H2Trees.noboundingsphereupdate,
            ),
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
        tree = H2Trees.buildtree(
            vertices(ms[i]);
            builder=H2Trees.KMeansTreeBuilder(;
                numberofclusters=10,
                minvalues=100,
                splitterkwargs=(; n_threads=1, rng=Random.MersenneTwister(1234)),
                updateradii=H2Trees.unsafemaxradiusboundingsphere,
            ),
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

    # Direct constructor with an all-defaults builder: deterministic (seeded `rng` default),
    # not just reachable via `buildtree`.
    default1 = H2Trees.KMeansTree(vertices(ms[1]))
    default2 = H2Trees.KMeansTree(vertices(ms[1]))
    @test default1 isa H2Trees.BoundingBallTree
    @test H2Trees.numberofnodes(default1) == H2Trees.numberofnodes(default2)

    X = raviartthomas(ms[1])
    _spacebuilder() = H2Trees.KMeansTreeBuilder(;
        numberofclusters=4,
        splitterkwargs=(; n_threads=1, rng=Random.MersenneTwister(1234)),
    )
    viaspace = H2Trees.KMeansTree(X; builder=_spacebuilder())
    viabuildtree = H2Trees.buildtree(X; builder=_spacebuilder())
    viapoints = H2Trees.buildtree(BEAST.positions(X); builder=_spacebuilder())
    @test viaspace isa H2Trees.BoundingBallTree
    @test H2Trees.numberofnodes(viaspace) == H2Trees.numberofnodes(viabuildtree)
    @test H2Trees.numberofnodes(viabuildtree) == H2Trees.numberofnodes(viapoints)
    @test_throws ArgumentError H2Trees.buildtree(
        X; graphweights=(nothing, nothing), builder=_spacebuilder()
    )

    @testset "balanceleaves! on realistic KMeans geometries" begin
        sawunbalanced = false
        for (i, m) in enumerate(ms)
            points = vertices(m)
            tree = H2Trees.buildtree(
                points;
                builder=H2Trees.KMeansTreeBuilder(;
                    numberofclusters=4,
                    minvalues=40,
                    splitterkwargs=(;
                        n_threads=1, rng=Random.MersenneTwister(20260803 + i)
                    ),
                ),
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
                for point in H2Trees.values(tree, leaf)
                    @test H2Trees.isin(tree, leaf, points[point])
                end
            end
        end
        @test sawunbalanced
    end

    @testset "empty k-means clusters don't crash construction" begin
        # `kmeans(...)` can converge with fewer populated clusters than requested (confirmed
        # here: with only 3 near-duplicate point locations and numberofclusters=10, it reliably
        # returns 10 centers but only ever assigns points to 3 of them, for every seed tried).
        # `kmeanswrapper` used to size `partitions` from `maximum(kresult.assignments)` while
        # `centers` came from `size(kresult.centers, 2)` -- whenever a requested cluster ends up
        # with zero points, those two counts diverge and indexing `partitions[i]` for the
        # unassigned cluster throws BoundsError. Fixed by sizing both from the cluster count and
        # filtering out empty partitions (mirroring `MetisTree.jl`'s `metiswrapper`, which faces
        # the same "a requested partition ended up empty" case).
        locations = [SVector(0.0, 0.0), SVector(0.001, 0.0), SVector(0.0, 0.001)]
        for seed in (142, 143, 144, 150, 175, 200)
            rng = Random.MersenneTwister(seed)
            points = [locations[rand(rng, 1:3)] for _ in 1:30]

            tree = H2Trees.buildtree(
                points;
                builder=H2Trees.KMeansTreeBuilder(;
                    numberofclusters=10,
                    minvalues=15,
                    splitterkwargs=(; n_threads=1, rng=Random.MersenneTwister(seed)),
                ),
            )

            @test sort!(reduce(vcat, H2Trees.values.(Ref(tree), H2Trees.leaves(tree)))) ==
                1:length(points)
            for leaf in H2Trees.leaves(tree)
                @test !isempty(H2Trees.values(tree, leaf))
            end
        end
    end
end

end # module TestH2ParallelKMeansTrees
