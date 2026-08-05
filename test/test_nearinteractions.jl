module TestNearInteractions

using Test
using CompScienceMeshes
using StaticArrays
using H2Trees
using BEAST
using SparseArrays

@testset "Near interaction argument validation" begin
    @test_throws(
        ArgumentError("η must be greater than or equal to 1, got 0.5"),
        H2Trees.isnearradius(SVector(0.0, 0.0), SVector(1.0, 1.0), 1.0, 1.0; η=0.5)
    )

    points = [SVector(0.0, 0.0), SVector(1.0, 0.0)]
    tree = buildtree(points; builder=TwoNTreeBuilder(; minhalfsize=0.25, minvalues=1))
    @test_throws(
        ArgumentError(
            "extractselfvalues is only supported for single-tree nearinteractions"
        ),
        H2Trees.nearinteractions(tree, tree; extractselfvalues=true)
    )
end

@testset "Galerkin Nearinteractions" begin
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

    for m in ms
        X = raviartthomas(m)

        minhalfsize = λ / 9
        tree = buildtree(
            X; builder=TwoNTreeBuilder(; minhalfsize=minhalfsize, minvalues=10)
        )

        blocktree = H2Trees.BlockTree(tree, tree)

        values, nearvalues = H2Trees.nearinteractions(tree;)

        selfvalues2, values2, nearvalues2 = H2Trees.nearinteractions(
            tree; extractselfvalues=true
        )

        testvalues, trialvalues = H2Trees.nearinteractions(blocktree)

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
                        push!(I, nv)
                        push!(J, v)
                    end
                end
            end
        end
        Atest = sparse(I, J, ones(length(I)), numfunctions(X), numfunctions(X))
        Atest.nzval .= 1.0

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
            tree = buildtree(
                X,
                Y;
                builder=BlockTreeBuilder(;
                    test=TwoNTreeBuilder(; minhalfsize=minhalfsize, minvalues=10),
                    trial=TwoNTreeBuilder(; minhalfsize=minhalfsize, minvalues=10),
                ),
            )
            testtree = H2Trees.testtree(tree)
            trialtree = H2Trees.trialtree(tree)

            values, nearvalues = H2Trees.nearinteractions(tree)

            @test length(values) == length(nearvalues)

            isempty(values) && continue
            isempty(nearvalues) && continue

            @show sum(length, values) / length(values)
            @show sum(length, nearvalues) / length(nearvalues)

            I = Int[]
            J = Int[]

            for leaf in H2Trees.leaves(testtree)
                leafvalues = H2Trees.values(testtree, leaf)
                _nearvalues = H2Trees.nearnodevalues(trialtree, testtree, leaf)
                for v in leafvalues
                    for nv in _nearvalues
                        push!(I, v)
                        push!(J, nv)
                    end
                end
                _farvalues = H2Trees.farnodevalues(trialtree, testtree, leaf)
                @test length(_nearvalues) + length(_farvalues) == numfunctions(Y)
            end

            Atest = sparse(I, J, ones(length(I)), numfunctions(X), numfunctions(Y))
            Atest.nzval .= 1.0
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

end # module TestNearInteractions
