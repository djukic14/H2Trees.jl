using Test
using CompScienceMeshes
using StaticArrays
using H2Trees
using BEAST

@testset "Empty tree" begin
    tree = TwoNTree(SVector(0.0, 0.0, 0.0), 1.0)

    @test H2Trees.root(tree) == 1
    @test H2Trees.center(tree, 1) == SVector(0.0, 0.0, 0.0)
    @test H2Trees.halfsize(tree, 1) == 1.0
    @test H2Trees.level(tree, 1) == 1

    @test H2Trees.LevelIterator(tree, 1) == Int[]
    @test H2Trees.treetrait(tree) == H2Trees.isTwoNTree()
end

@testset "Filled Tree" begin
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere4.in")
    )

    points = vertices(m)

    root = 2
    minlevel = 2

    tree = TwoNTree(points, 0.0; root=root, minlevel=minlevel, minvalues=10)
    valuesatnodes = H2Trees.valuesatnodes(tree)
    @test length(valuesatnodes) == length(points)
    for (functionid, value) in enumerate(valuesatnodes)
        @test length(value) == 1
        @test functionid in H2Trees.values(tree, value[1])
    end
    nodesatvalues = H2Trees.nodesatvalues(tree)
    for (key, value) in nodesatvalues
        @test length(key) == 1
        key = key[1]
        @test sort(value) == sort(H2Trees.values(tree, key))
    end

    maximumlevel = H2Trees.levels(tree)[end] - minlevel + 1

    leaflevels = sort(unique(H2Trees.level.(Ref(tree), H2Trees.leaves(tree))))[2:end]

    nodes = Int[]
    for i in H2Trees.DepthFirstIterator(tree, root)
        @test any(i .== H2Trees.LevelIterator(tree, H2Trees.level(tree, i)))

        for node in H2Trees.SameLevelIterator(tree, i)
            @test H2Trees.level(tree, node) == H2Trees.level(tree, i)
        end

        for node in H2Trees.TranslatingNodesIterator(tree, i)
            @test H2Trees.level(tree, node) == H2Trees.level(tree, i)
            @test !H2Trees.isnear(tree, node, i)
        end

        valuesonlevel = Int[]

        for node in H2Trees.TranslatingNodesIterator(tree, i)
            append!(valuesonlevel, H2Trees.values(tree, node))
        end

        for node in H2Trees.NotTranslatingNodesIterator(tree, i)
            append!(valuesonlevel, H2Trees.values(tree, node))
        end

        if !(H2Trees.level(tree, i) in leaflevels)
            @test sort!(valuesonlevel) == Array(1:length(points))
        end
        push!(nodes, i)

        sector = H2Trees.sector(tree, i)
        oppositesector = H2Trees.oppositesector(tree, i)

        if sector == 0
            @test oppositesector == 7
        elseif sector == 1
            @test oppositesector == 6
        elseif sector == 2
            @test oppositesector == 5
        elseif sector == 3
            @test oppositesector == 4
        elseif sector == 4
            @test oppositesector == 3
        elseif sector == 5
            @test oppositesector == 2
        elseif sector == 6
            @test oppositesector == 1
        elseif sector == 7
            @test oppositesector == 0
        end

        @test H2Trees.levelindex(tree, i) == H2Trees.level(tree, i) - minlevel + 1
        i == root && continue

        pminuschild = H2Trees.parentcenterminuschildcenter(tree, i)

        @test pminuschild ≈
            H2Trees.center(tree, H2Trees.parent(tree, i)) - H2Trees.center(tree, i)

        for cornerid in 1:8
            corner = H2Trees.cornerpoints(tree, i, cornerid)

            if cornerid == 1
                @test corner ≈ H2Trees.center(tree, i) .- H2Trees.halfsize(tree, i)

            elseif cornerid == 2
                @test corner ≈
                    H2Trees.center(tree, i) +
                      H2Trees.halfsize(tree, i) .* SVector(-1, -1, 1)

            elseif cornerid == 3
                @test corner ≈
                    H2Trees.center(tree, i) +
                      H2Trees.halfsize(tree, i) .* SVector(-1, 1, -1)

            elseif cornerid == 4
                @test corner ≈
                    H2Trees.center(tree, i) +
                      H2Trees.halfsize(tree, i) .* SVector(-1, 1, 1)

            elseif cornerid == 5
                @test corner ≈
                    H2Trees.center(tree, i) +
                      H2Trees.halfsize(tree, i) .* SVector(1, -1, -1)
            elseif cornerid == 6
                @test corner ≈
                    H2Trees.center(tree, i) +
                      H2Trees.halfsize(tree, i) .* SVector(1, -1, 1)

            elseif cornerid == 7
                @test corner ≈
                    H2Trees.center(tree, i) +
                      H2Trees.halfsize(tree, i) .* SVector(1, 1, -1)
            else
                cornerid == 8
                @test corner ≈
                    H2Trees.center(tree, i) +
                      H2Trees.halfsize(tree, i) .* SVector(1, 1, 1)
            end
        end
    end

    @test sort!(nodes) == Array(root:(length(tree.nodes) + root - 1))

    for level in H2Trees.levels(tree)
        level in leaflevels && continue
        values = Int[]
        for node in H2Trees.LevelIterator(tree, level)
            append!(values, H2Trees.values(tree, node))
        end

        @test sort(values) == Array(1:length(points))
    end

    leaves = H2Trees.leaves(tree)

    for leaf in leaves
        @test H2Trees.isleaf(tree, leaf)
        @test tree(leaf).firstchild == 0
    end

    tree2 = TwoNTree(SVector(0.0, 0.0, 0.0), 0.1; root=root, minlevel=minlevel)
    @test H2Trees.root(tree2) == root
    @test H2Trees.center(tree2, root) == SVector(0.0, 0.0, 0.0)
    @test H2Trees.halfsize(tree2, root) == 0.1
    @test H2Trees.level(tree2, root) == minlevel

    @test H2Trees.treetrait(tree) == H2Trees.isTwoNTree()
    @test H2Trees.treetrait(tree2) == H2Trees.isTwoNTree()

    @test H2Trees.minhalfsize(tree) ≈ H2Trees.halfsize(tree) * (1 / 2)^(maximumlevel - 1)

    for value in eachindex(points)
        @test H2Trees.findleafnode(tree, value) ∈ leaves
        @test value ∈ H2Trees.values(tree, H2Trees.findleafnode(tree, value))
    end

    @test H2Trees.findleafnode(tree, length(points) + 1) == 0
    @test H2Trees.findleafnode(tree, -1) == 0

    leafclusters = H2Trees.leafclusters(tree)

    for (i, leaf) in enumerate(H2Trees.leaves(tree))
        @test leafclusters[i] == H2Trees.values(tree, leaf)
    end
end

@testset "Float32" begin
    tree = TwoNTree(SVector(0.0f0, 0.0f0, 0.0f0), 1.0f0)

    @test H2Trees.root(tree) == 1
    @test H2Trees.center(tree, 1) == SVector(0.0f0, 0.0f0, 0.0f0)
    @test H2Trees.halfsize(tree, 1) == 1.0f0
    @test H2Trees.level(tree, 1) == 1

    @test H2Trees.treetrait(tree) == H2Trees.isTwoNTree()
end

@testset "Simple Blocktree" begin
    λ = 1.0
    minhalfsize = λ / 10

    mx = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere5.in")
    )

    my = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter7.in")
    )

    X = raviartthomas(mx)
    Y = raviartthomas(my)

    tree = TwoNTree(X, Y, minhalfsize; minvaluestest=10, minvaluestrial=3)

    for tree in [H2Trees.testtree(tree), H2Trees.trialtree(tree)]
        valuesatnodes = H2Trees.valuesatnodes(tree)
        for (functionid, value) in enumerate(valuesatnodes)
            @test length(value) == 1
            @test functionid in H2Trees.values(tree, value[1])
        end
        nodesatvalues = H2Trees.nodesatvalues(tree)
        for (key, value) in nodesatvalues
            @test length(key) == 1
            key = key[1]
            @test sort(value) == sort(H2Trees.values(tree, key))
        end
    end

    @test eltype(tree) == SVector{3,Float64}
    @test eltype(H2Trees.testtree(tree)) == SVector{3,Float64}

    @test H2Trees.treewithmorelevels(tree) == H2Trees.trialtree(tree)

    tree2 = TwoNTree(Y, X, minhalfsize)
    @test H2Trees.treewithmorelevels(tree2) == H2Trees.trialtree(tree2)

    @test H2Trees.minhalfsize(H2Trees.trialtree(tree)) == minhalfsize

    for level in H2Trees.levels(H2Trees.trialtree(tree))
        halfsize = H2Trees.halfsize(
            H2Trees.trialtree(tree),
            H2Trees.LevelIterator(H2Trees.trialtree(tree), level)[begin],
        )

        for node in H2Trees.LevelIterator(H2Trees.trialtree(tree), level)
            @test H2Trees.level(H2Trees.trialtree(tree), node) == level
            @test H2Trees.halfsize(H2Trees.trialtree(tree), node) == halfsize

            for testnode in H2Trees.LevelIterator(H2Trees.testtree(tree), level)
                @test H2Trees.level(H2Trees.trialtree(tree), node) ==
                    H2Trees.level(H2Trees.testtree(tree), testnode)
                @test H2Trees.halfsize(H2Trees.trialtree(tree), node) ==
                    H2Trees.halfsize(H2Trees.testtree(tree), testnode)
            end
        end

        level == 1 && continue

        level - 1 ∉ H2Trees.levels(H2Trees.trialtree(tree)) && continue

        halfsizeabove = H2Trees.halfsize(
            H2Trees.trialtree(tree),
            H2Trees.LevelIterator(H2Trees.trialtree(tree), level - 1)[begin],
        )

        @test halfsizeabove == 2 * halfsize
    end

    @test H2Trees.treetrait(tree) == H2Trees.isBlockTree()

    display(tree)
    display(H2Trees.testtree(tree))
    display(H2Trees.trialtree(tree))
end
