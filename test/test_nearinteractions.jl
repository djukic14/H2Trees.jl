using Test
using CompScienceMeshes
using StaticArrays
using H2Trees
using BEAST
using SparseArrays

@testset "Galerkin Nearinteractions" begin
    λ = 1.0

    ms = [
        CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere3.in")
        ),
        # CompScienceMeshes.readmesh(
        #     joinpath(pkgdir(H2Trees), "test", "assets", "in", "cuboid.in")
        # ),
        # CompScienceMeshes.readmesh(
        #     joinpath(pkgdir(H2Trees), "test", "assets", "in", "twospheres2.in")
        # ),
        CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "multiplerects.in")
        ),
    ]

    for m in ms
        X = raviartthomas(m)

        minhalfsize = λ / 9
        tree = TwoNTree(X, minhalfsize)

        blocktree = H2Trees.BlockTree(tree, tree)

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

        for i in eachindex(values)
            for v in values[i]
                for nv in nearvalues[i]
                    A[v, nv] += 1.0
                end
            end
        end

        @test all(nonzeros(A) .== 2)

        for i in eachindex(values2)
            for v in values2[i]
                for nv in nearvalues2[i]
                    A2[v, nv] += 1.0
                end
            end
        end

        for i in eachindex(selfvalues2)
            for v in selfvalues2[i]
                for nv in selfvalues2[i]
                    A2[v, nv] += 1.0
                end
            end
        end

        @test all(nonzeros(A2) .== 2)

        for i in eachindex(testvalues)
            for v in testvalues[i]
                for nv in trialvalues[i]
                    Ablock[v, nv] += 1.0
                end
            end
        end

        @test all(nonzeros(Ablock) .== 2)
    end
end

@testset "Petrov-Galerkin Nearinteractions" begin
    λ = 1.0

    ms = [
        CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere3.in")
        ),
        CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "cuboid.in")
        ),
        CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "twospheres2.in")
        ),
        CompScienceMeshes.readmesh(
            joinpath(pkgdir(H2Trees), "test", "assets", "in", "multiplerects.in")
        ),
    ]

    for mx in ms
        for my in ms
            X = raviartthomas(mx)
            Y = raviartthomas(my)

            minhalfsize = λ / 9
            tree = TwoNTree(X, Y, minhalfsize)
            testtree = H2Trees.testtree(tree)
            trialtree = H2Trees.trialtree(tree)

            values, nearvalues = H2Trees.nearinteractions(tree; size=30)

            @test length(values) == length(nearvalues)

            isempty(values) && continue
            isempty(nearvalues) && continue

            @show sum(length, values) / length(values)
            @show sum(length, nearvalues) / length(nearvalues)

            I = Int[]
            J = Int[]

            for leaf in H2Trees.leaves(testtree)
                leafvalues = H2Trees.values(testtree, leaf)
                for node in H2Trees.NearNodeIterator(trialtree, testtree, leaf)
                    nearnodevalues = H2Trees.values(trialtree, node)
                    for v in leafvalues
                        for nv in nearnodevalues
                            push!(I, v)
                            push!(J, nv)
                        end
                    end
                end
            end
            Atest = sparse(I, J, ones(length(I)), numfunctions(X), numfunctions(Y))

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

            A = sparse(I, J, ones(length(I)), numfunctions(X), numfunctions(Y))

            @test maximum(abs, A - Atest) == 0

            for i in eachindex(values)
                for v in values[i]
                    for nv in nearvalues[i]
                        A[v, nv] += 1.0
                    end
                end
            end

            @test all(nonzeros(A) .== 2)
        end
    end
end
