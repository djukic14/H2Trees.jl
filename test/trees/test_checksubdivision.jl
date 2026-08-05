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
                            tree = buildtree(
                                X;
                                builder=TwoNTreeBuilder(;
                                    minlevel=minlevel,
                                    root=root,
                                    minvalues=minvalues,
                                    minhalfsize=minhalfsize,
                                    protrusion=ProtrusionCheck(;
                                        max=maxprotrusion, compute=protrusion
                                    ),
                                ),
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

@testset "Block TwoNTree Check Subdivision" begin
    meshes = [
        "cuboid",
        # "cuboid2",
        # "cuboid3",
        # "cuboid4",
        "multiplerects",
        # "multiplerects2",
        # "multiplerects3",
        # "sphere",
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

    for (ix, mx) in enumerate(ms)
        X = raviartthomas(mx)
        testprotrusion = H2Trees.BEASTProtrusionFunctor(X)
        for (iy, my) in enumerate(ms)
            Y = raviartthomas(my)
            trialprotrusion = H2Trees.BEASTProtrusionFunctor(Y)
            for testmaxprotrusion in [0.1, 0.5]
                for trialmaxprotrusion in [0.1, 0.5]
                    for testminvalues in [0, 100]
                        for trialminvalues in [0, 100]
                            for minhalfsize in [0.0, 0.1]
                                println(
                                    "Testing with meshes $ix and $iy, testmaxprotrusion=$testmaxprotrusion, trialmaxprotrusion=$trialmaxprotrusion, testminvalues=$testminvalues, trialminvalues=$trialminvalues, minhalfsize=$minhalfsize",
                                )
                                tree = buildtree(
                                    X,
                                    Y;
                                    builder=BlockTreeBuilder(;
                                        test=TwoNTreeBuilder(;
                                            minhalfsize=minhalfsize,
                                            minvalues=testminvalues,
                                            protrusion=ProtrusionCheck(;
                                                max=testmaxprotrusion,
                                                compute=testprotrusion,
                                            ),
                                        ),
                                        trial=TwoNTreeBuilder(;
                                            minhalfsize=minhalfsize,
                                            minvalues=trialminvalues,
                                            protrusion=ProtrusionCheck(;
                                                max=trialmaxprotrusion,
                                                compute=trialprotrusion,
                                            ),
                                        ),
                                    ),
                                )
                                testtree = H2Trees.testtree(tree)
                                trialtree = H2Trees.trialtree(tree)

                                _testmaxprotrusion = H2Trees.maxprotrusion(
                                    testtree; computeprotrusion=testprotrusion
                                )
                                @test maximum(_testmaxprotrusion) <= testmaxprotrusion

                                _trialmaxprotrusion = H2Trees.maxprotrusion(
                                    trialtree; computeprotrusion=trialprotrusion
                                )
                                @test maximum(_trialmaxprotrusion) <= trialmaxprotrusion

                                @test H2Trees.minhalfsize(testtree) >= minhalfsize
                                @test H2Trees.minhalfsize(trialtree) >= minhalfsize

                                for leaf in H2Trees.leaves(testtree)
                                    leaf == H2Trees.root(testtree) && continue
                                    @test length(
                                        H2Trees.values(
                                            testtree, H2Trees.parent(testtree, leaf)
                                        ),
                                    ) >= testminvalues
                                end

                                for leaf in H2Trees.leaves(trialtree)
                                    leaf == H2Trees.root(trialtree) && continue
                                    @test length(
                                        H2Trees.values(
                                            trialtree, H2Trees.parent(trialtree, leaf)
                                        ),
                                    ) >= trialminvalues
                                end

                                tesths = H2Trees.halfsizes(testtree)
                                trialhs = H2Trees.halfsizes(trialtree)

                                for (i, level) in enumerate(H2Trees.levels(testtree))
                                    !(level in H2Trees.levels(trialtree)) && continue
                                    j = findfirst(==(level), H2Trees.levels(trialtree))

                                    @test tesths[i] ≈ trialhs[j]
                                end

                                for (i, level) in enumerate(H2Trees.levels(trialtree))
                                    !(level in H2Trees.levels(testtree)) && continue
                                    j = findfirst(==(level), H2Trees.levels(testtree))

                                    @test trialhs[i] ≈ tesths[j]
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

@testset "BEAST spaces resolve automatic protrusion defaults" begin
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "cuboid.in")
    )
    X = raviartthomas(m)

    # Plain point collections keep the builder default as a no-op protrusion policy.
    pointsbuilder = H2Trees._resolve_builder_protrusion(
        TwoNTreeBuilder(), BEAST.positions(X)
    )
    @test pointsbuilder.protrusion == NoProtrusionCheck()

    # BEAST spaces turn the automatic policy into the old space-constructor default: element
    # protrusions with a conservative 0.25 threshold.
    spacebuilder = H2Trees._resolve_builder_protrusion(TwoNTreeBuilder(; minvalues=120), X)
    @test spacebuilder.protrusion.compute isa H2Trees.BEASTProtrusionFunctor
    @test spacebuilder.protrusion.max == 0.25
    @test spacebuilder.minvalues == 120

    # Explicit max overrides keep the requested threshold while still using the BEAST functor.
    override = TwoNTreeBuilder(; minvalues=120, protrusion=ProtrusionCheck(; max=0.5))
    adjusted = H2Trees._resolve_builder_protrusion(override, X)
    @test adjusted.protrusion.compute isa H2Trees.BEASTProtrusionFunctor
    @test adjusted.protrusion.max == 0.5

    # The direct constructor and the canonical builder path are the same API surface.
    directtree = TwoNTree(X; builder=TwoNTreeBuilder(; minvalues=120))
    autotree = buildtree(X; builder=TwoNTreeBuilder(; minvalues=120))
    explicittree = buildtree(
        X;
        builder=TwoNTreeBuilder(;
            minvalues=120,
            protrusion=ProtrusionCheck(;
                max=0.25, compute=H2Trees.BEASTProtrusionFunctor(X)
            ),
        ),
    )
    @test H2Trees.leaves(directtree) == H2Trees.leaves(autotree)
    @test H2Trees.leaves(autotree) == H2Trees.leaves(explicittree)

    # BlockTree(testspace, trialspace; builder) is the same API surface as buildtree, mirroring the
    # single-space TwoNTree/buildtree parity checked above -- including that AutoProtrusionCheck
    # resolves to the BEAST default on the direct-constructor path too, not only through buildtree.
    blockbuilder = BlockTreeBuilder(;
        test=TwoNTreeBuilder(; minvalues=120), trial=TwoNTreeBuilder(; minvalues=120)
    )
    explicitblockbuilder = BlockTreeBuilder(;
        test=TwoNTreeBuilder(;
            minvalues=120,
            protrusion=ProtrusionCheck(;
                max=0.25, compute=H2Trees.BEASTProtrusionFunctor(X)
            ),
        ),
        trial=TwoNTreeBuilder(;
            minvalues=120,
            protrusion=ProtrusionCheck(;
                max=0.25, compute=H2Trees.BEASTProtrusionFunctor(X)
            ),
        ),
    )
    directblock = BlockTree(X, X; builder=blockbuilder)
    autoblock = buildtree(X, X; builder=blockbuilder)
    explicitblock = buildtree(X, X; builder=explicitblockbuilder)
    @test H2Trees.leaves(H2Trees.testtree(directblock)) ==
        H2Trees.leaves(H2Trees.testtree(autoblock)) ==
        H2Trees.leaves(H2Trees.testtree(explicitblock))
    @test H2Trees.leaves(H2Trees.trialtree(directblock)) ==
        H2Trees.leaves(H2Trees.trialtree(autoblock)) ==
        H2Trees.leaves(H2Trees.trialtree(explicitblock))

    @test_throws ArgumentError buildtree(
        X; graphweights=(nothing, nothing), builder=TwoNTreeBuilder()
    )

    # Explicit opt-out and non-`TwoNTreeBuilder` builders are left untouched.
    noprot = H2Trees._resolve_builder_protrusion(
        TwoNTreeBuilder(; minhalfsize=0.1, protrusion=NoProtrusionCheck()), X
    )
    @test noprot.protrusion == NoProtrusionCheck()
    @test H2Trees._resolve_builder_protrusion(KMeansTreeBuilder(), X) isa KMeansTreeBuilder
end
