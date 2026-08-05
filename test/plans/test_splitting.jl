using Test
using BEAST, CompScienceMeshes
using StaticArrays
using LinearAlgebra
using H2Trees

@testset "Spit Aggregate Translate Plan" begin
    λ = 1.0
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere6.in")
    )
    X = raviartthomas(m)
    tree = buildtree(
        buildtree(
            X;
            builder=TwoNTreeBuilder(;
                minhalfsize=λ / 10, minvalues=0, protrusion=NoProtrusionCheck()
            ),
        );
        builder=SimpleHybridTreeBuilder(; hybridhalfsize=λ / 5),
    )

    TFIterator = H2Trees.TranslatingNodesIterator(; isnear=H2Trees.isnear())(tree)
    aggregatetranslateplan = H2Trees.AggregateTranslatePlan(tree, TFIterator)
    upperplan, lowerplan = H2Trees.splitplan(tree, aggregatetranslateplan)

    @test H2Trees.levels(upperplan) == H2Trees.hybridlevel(tree):-1:3
    @test H2Trees.levels(lowerplan) ==
        H2Trees.levels(tree)[end]:-1:(H2Trees.hybridlevel(tree) + 1)

    @test H2Trees.tree(upperplan) === tree
    @test H2Trees.tree(lowerplan) === tree
    for level in H2Trees.levels(upperplan)
        @test H2Trees.receivingnodes(upperplan, level) ===
            upperplan.receivingnodes_by_level[H2Trees.leveltolevelid(upperplan, level)]
    end
    for level in H2Trees.levels(lowerplan)
        @test H2Trees.receivingnodes(lowerplan, level) ===
            lowerplan.receivingnodes_by_level[H2Trees.leveltolevelid(lowerplan, level)]
    end

    for nodes in H2Trees.nodes(upperplan)
        for node in nodes
            @test H2Trees.isuppertreenode(tree, node)
        end
    end

    for nodes in H2Trees.nodes(lowerplan)
        for node in nodes
            @test H2Trees.islowertreenode(tree, node)
        end
    end

    for receivingnodes in H2Trees.receivingnodes(upperplan)
        for (node, receivingnodes) in receivingnodes
            @test H2Trees.isuppertreenode(tree, node)
            for node in receivingnodes
                @test H2Trees.isuppertreenode(tree, node)
            end
        end
    end

    for receivingnodes in H2Trees.receivingnodes(lowerplan)
        for (node, receivingnodes) in receivingnodes
            @test H2Trees.islowertreenode(tree, node)
            for node in receivingnodes
                @test H2Trees.islowertreenode(tree, node)
            end
        end
    end
end

@testset "Disaggregate Translate Plan" begin
    λ = 1.0
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere6.in")
    )
    X = raviartthomas(m)
    tree = buildtree(
        buildtree(
            X;
            builder=TwoNTreeBuilder(;
                minhalfsize=λ / 10, minvalues=0, protrusion=NoProtrusionCheck()
            ),
        );
        builder=SimpleHybridTreeBuilder(; hybridhalfsize=λ / 5),
    )

    disaggregationplan = H2Trees.DisaggregateTranslatePlan(
        tree, H2Trees.TranslatingNodesIterator(; isnear=H2Trees.isnear())(tree)
    )

    upperplan, lowerplan = H2Trees.splitplan(tree, disaggregationplan)

    @test upperplan.isdisaggregationnode == disaggregationplan.isdisaggregationnode
    @test lowerplan.isdisaggregationnode == disaggregationplan.isdisaggregationnode
    for level in H2Trees.levels(upperplan)
        @test H2Trees.receivingnodes(upperplan, level) ===
            upperplan.receivingnodes_by_level[H2Trees.leveltolevelid(upperplan, level)]
    end
    for level in H2Trees.levels(lowerplan)
        @test H2Trees.receivingnodes(lowerplan, level) ===
            lowerplan.receivingnodes_by_level[H2Trees.leveltolevelid(lowerplan, level)]
    end

    @test H2Trees.levels(upperplan) == 3:H2Trees.hybridlevel(tree)
    @test H2Trees.levels(lowerplan) ==
        (H2Trees.hybridlevel(tree) + 1):H2Trees.levels(tree)[end]

    for nodes in H2Trees.nodes(upperplan)
        for node in nodes
            @test H2Trees.isuppertreenode(tree, node)
        end
    end

    for nodes in H2Trees.nodes(lowerplan)
        for node in nodes
            @test H2Trees.islowertreenode(tree, node)
        end
    end

    for receivingnodes in H2Trees.translatingnodes(upperplan)
        for (node, receivingnodes) in receivingnodes
            @test H2Trees.isuppertreenode(tree, node)
            for node in receivingnodes
                @test H2Trees.isuppertreenode(tree, node)
            end
        end
    end

    for receivingnodes in H2Trees.translatingnodes(lowerplan)
        for (node, receivingnodes) in receivingnodes
            @test H2Trees.islowertreenode(tree, node)
            for node in receivingnodes
                @test H2Trees.islowertreenode(tree, node)
            end
        end
    end
end
