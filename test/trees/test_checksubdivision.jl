using Test
using CompScienceMeshes, BEAST
using H2Trees

@testset "TwoNTree Check Subdivision" begin
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
    # end

    m = ms[1]
    X = raviartthomas(m)
    protrusion = H2Trees.BEASTProtrusionFunctor(X)
    for m in ms
        for root in 1:2
            for minlevel in 1:2
                for maxprotrusion in [0.1, 0.5, 1.0]
                    for minvalues in [0, 10, 100]
                        for minhalfsize in [0.0, 0.1]
                            tree = TwoNTree(
                                X,
                                minhalfsize;
                                minlevel=minlevel,
                                root=root,
                                minvalues=minvalues,
                                maxprotrusion=maxprotrusion,
                                computeprotrusion=protrusion,
                            )

                            _maxprotrusion = H2Trees.maxprotrusion(
                                tree; computeprotrusion=protrusion
                            )
                            @test maximum(_maxprotrusion) <= maxprotrusion
                            @test H2Trees.minhalfsize(tree) >= minhalfsize
                            for leaf in H2Trees.leaves(tree)
                                @test length(
                                    H2Trees.values(tree, H2Trees.parent(tree, leaf))
                                ) >= minvalues
                            end
                        end
                    end
                end
            end
        end
    end
end
