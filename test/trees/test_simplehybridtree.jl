using Test
using CompScienceMeshes
using StaticArrays
using H2Trees

@testset "Simple Hybrid Tree" begin
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere6.in")
    )

    points = vertices(m)

    tree = TwoNTree(points, 0.1; root=2, minlevel=2, minvalues=10)

    @test_throws ErrorException hybridtree = SimpleHybridTree(tree; hybridhalfsize=0.1)

    @test_throws ErrorException hybridtree = SimpleHybridTree(tree; hybridhalfsize=0.01)

    hybridtree = SimpleHybridTree(tree; hybridhalfsize=0.2)
    @test H2Trees.hybridlevel(hybridtree) == 4
    @test H2Trees.uniquepointstreetrait(hybridtree) == H2Trees.uniquepointstreetrait(tree)

    @test_nowarn println(hybridtree)
    @test_nowarn display(hybridtree)
    @test_nowarn show(hybridtree)

    @test H2Trees.numberoflevels(hybridtree) == H2Trees.numberoflevels(tree)
    @test H2Trees.numberofnodes(hybridtree) == H2Trees.numberofnodes(tree)
    @test eltype(hybridtree) == eltype(tree)
    @test H2Trees.treetrait(hybridtree) == H2Trees.treetrait(tree)

    hybridlevel = 4
    @test H2Trees.hybridlevel(hybridtree) == hybridlevel
    for node in H2Trees.DepthFirstIterator(hybridtree)
        if H2Trees.level(hybridtree, node) <= hybridlevel
            @test H2Trees.isuppertreenode(hybridtree, node)
        else
            @test H2Trees.islowertreenode(hybridtree, node)
        end
    end

    for node in H2Trees.DepthFirstIterator(tree)
        @test H2Trees.samelevelnodes(hybridtree, node) == H2Trees.samelevelnodes(tree, node)
        @test H2Trees.values(hybridtree, node) == H2Trees.values(tree, node)
        @test H2Trees.sector(hybridtree, node) == H2Trees.sector(tree, node)
        @test H2Trees.parent(hybridtree, node) == H2Trees.parent(tree, node)
        @test H2Trees.nextsibling(hybridtree, node) == H2Trees.nextsibling(tree, node)
        @test H2Trees.parentcenterminuschildcenter(hybridtree, node) ==
            H2Trees.parentcenterminuschildcenter(tree, node)
        @test H2Trees.oppositesector(hybridtree, node) == H2Trees.oppositesector(tree, node)
    end
end
