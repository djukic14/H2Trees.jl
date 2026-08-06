using Test
using CompScienceMeshes
using StaticArrays
using H2Trees

@testset "_uniformseparationdepth rejects coincident points" begin
    # Two bit-identical points always land in the same octree sector no matter how far the
    # bisection recurses: `sectorcentersize`'s split only compares `pt .> ct`, which never differs
    # between them, so `center` converges toward the shared point and `halfsize` halves toward
    # 0.0 forever without ever separating them. Before the fix this ran until the process stack
    # overflowed instead of returning or erroring.
    duplicate = SVector(0.0, 0.0, 0.0)
    @test_throws ArgumentError H2Trees._uniformseparationdepth(
        [duplicate, duplicate], [1, 2], duplicate, 1.0
    )

    # Three-way and off-origin coincidences hit the same code path.
    other = SVector(1.5, -2.0, 0.25)
    @test_throws ArgumentError H2Trees._uniformseparationdepth(
        [other, other, other], [1, 2, 3], other, 3.0
    )

    # The general adaptive `bulkbuildtwontree` path is the only construction path that ever calls
    # `_uniformseparationdepth` (the `minvalues=0`+`NoProtrusionCheck` uniform-bulk path
    # deliberately never does, see its docstring) -- it must surface the same clear error rather
    # than hanging/crashing.
    @test_throws ArgumentError H2Trees.buildtree(
        [duplicate, duplicate]; builder=H2Trees.TwoNTreeBuilder(; minvalues=1)
    )

    # Distinct points -- even very close ones -- must still separate normally and are unaffected.
    close = SVector(1e-10, 0.0, 0.0)
    @test H2Trees._uniformseparationdepth([duplicate, close], [1, 2], duplicate, 1.0) > 0
end

@testset "locatepoint on a uniform-depth tree" begin
    # `locatepoint` walks from the root toward a point and stops once it reaches the requested
    # level -- it requires every branch to actually reach that level. The `minvalues=0`+
    # `NoProtrusionCheck` uniform-bulk path guarantees this (recursion is bounded purely by
    # `halfsize <= minhalfsize`, identical along every branch regardless of point distribution),
    # unlike the general adaptive path, where leaves can stop at different depths.
    meshes =
        ["cuboid", "multiplerects", "sphere", "spherewithcenter", "twospheres"] .* ".in"

    ms = [
        CompScienceMeshes.readmesh(joinpath(pkgdir(H2Trees), "test", "assets", "in", m)) for
        m in meshes
    ]

    for minlevel in 1:3
        for root in 1:2
            for (i, m) in enumerate(ms)
                points = vertices(m)
                _, rootsize = H2Trees.boundingbox(points)
                tree = H2Trees.buildtree(
                    points;
                    builder=H2Trees.TwoNTreeBuilder(;
                        minlevel=minlevel,
                        root=root,
                        minvalues=0,
                        minhalfsize=rootsize / 16,
                        protrusion=H2Trees.NoProtrusionCheck(),
                    ),
                )

                @test H2Trees.checkbalancedtree(tree)
                @test H2Trees.levels(tree)[1] == minlevel
                @test H2Trees.root(tree) == root
                @test length(H2Trees.levels(tree)) > 1

                for leaf in H2Trees.leaves(tree)
                    for val in H2Trees.values(tree, leaf)
                        @test H2Trees.isin(tree, leaf, points[val])
                    end
                end

                for level in H2Trees.levels(tree)
                    for point in points
                        node = H2Trees.locatepoint(tree, point, level)
                        @test H2Trees.isin(tree, node, point)
                        @test H2Trees.level(tree, node) == level
                    end
                end
                @test_throws ErrorException H2Trees.locatepoint(
                    tree, points[1], H2Trees.levels(tree)[end] + 1
                )
            end
        end
    end
end
