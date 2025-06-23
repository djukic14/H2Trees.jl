using CompScienceMeshes, BEAST
using H2Trees
using Test
using SparseArrays

@testset "QuadPointsTree" begin
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

    for m in ms
        X = raviartthomas(m)

        minhalfsize = λ / 10
        tree = QuadPointsTree(X, minhalfsize)

        @test H2Trees.numberofvalues(tree) == numfunctions(X)

        sort!(unique!(H2Trees.values(tree, H2Trees.root(tree)))) == 1:numfunctions(X)

        valuesatnodes = H2Trees.valuesatnodes(tree)
        for (functionid, nodes) in enumerate(valuesatnodes)
            @test issorted(nodes)
            for node in nodes
                @test functionid in H2Trees.values(tree, node)
            end
        end

        nodesatvalues = H2Trees.nodesatvalues(tree)

        for (nodes, values) in nodesatvalues
            @test issorted(nodes)
            for node in nodes
                for value in values
                    @test value in H2Trees.values(tree, node)
                end
            end

            for value in values
                @test valuesatnodes[value] == nodes
            end
        end
    end
end

@testset "QuadPointsTree Nearinteractions" begin
    λ = 1.0

    # ms = [
    #     meshsphere(2 * λ, λ / 10),
    #     meshcuboid(λ, 2 * λ, 3 * λ, λ / 10),
    #     meshtwospheres([0.0, 0, 0], λ, [3 * λ, 0, 0], λ / 2, λ / 10),
    #     meshmultiplerectangles(λ, λ / 10, 4),
    # ]

    ms = [
        CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere8.in")
        ),
        CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "cuboid4.in")
        ),
        CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "twospheres4.in")
        ),
        CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "multiplerects3.in")
        ),
    ]

    for m in ms
        X = raviartthomas(m)

        minhalfsize = λ / 9
        tree = QuadPointsTree(X, minhalfsize)

        blocktree = BlockTree(tree, tree)

        values, nearvalues = H2Trees.nearinteractions(tree; size=30)

        selfvalues2, values2, nearvalues2 = H2Trees.nearinteractions(
            tree; size=30, extractselfvalues=true
        )

        testvalues, trialvalues = H2Trees.nearinteractions(blocktree; size=30)

        @show sum(length, values) / length(values)
        @show sum(length, nearvalues) / length(nearvalues)

        @show sum(length, values2) / length(values2)
        @show sum(length, nearvalues2) / length(nearvalues2)

        @show sum(length, testvalues) / length(testvalues)
        @show sum(length, trialvalues) / length(trialvalues)

        I = Int[]
        J = Int[]

        for leaf in H2Trees.leaves(tree)
            leafvalues = H2Trees.values(tree, leaf)
            for node in H2Trees.NearNodeIterator(tree, leaf)
                nearnodevalues = H2Trees.values(tree, node)
                for v in leafvalues
                    for nv in nearnodevalues
                        push!(I, v)
                        push!(J, nv)
                    end
                end
            end
        end

        Atest = sparse(I, J, ones(length(I)), numfunctions(X), numfunctions(X))
        nonzeros(Atest) .= 1.0

        @test length(values) == length(nearvalues)

        I = Int[]
        J = Int[]

        for i in eachindex(values)
            for v in values[i]
                for nv in nearvalues[i]
                    push!(I, v)
                    push!(J, nv)
                end
            end
        end

        A = sparse(I, J, ones(length(I)), numfunctions(X), numfunctions(X))

        @test maximum(abs, A - Atest) == 0

        I = Int[]
        J = Int[]

        for i in eachindex(values2)
            for v in values2[i]
                for nv in nearvalues2[i]
                    push!(I, v)
                    push!(J, nv)
                end
            end
        end

        for i in eachindex(selfvalues2)
            for v in selfvalues2[i]
                for nv in selfvalues2[i]
                    push!(I, v)
                    push!(J, nv)
                end
            end
        end
        A2 = sparse(I, J, ones(length(I)), numfunctions(X), numfunctions(X))

        @test maximum(abs, A2 - Atest) == 0

        I = Int[]
        J = Int[]

        for i in eachindex(testvalues)
            for v in testvalues[i]
                for nv in trialvalues[i]
                    push!(I, v)
                    push!(J, nv)
                end
            end
        end

        Ablock = sparse(I, J, ones(length(I)), numfunctions(X), numfunctions(X))
        @test maximum(abs, Ablock - Atest) == 0
    end
end
