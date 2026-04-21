using Test
using BEAST, CompScienceMeshes
using Metis, Graphs
using H2Trees

@testset "Adjacency graphs" begin
    meshes =
        ["cuboid", "multiplerects", "sphere", "spherewithcenter", "twospheres"] .* ".in"

    ms = [
        CompScienceMeshes.readmesh(joinpath(pkgdir(H2Trees), "test", "assets", "in", m)) for
        m in meshes
    ]

    m = ms[1]
    for (i, m) in enumerate(ms)
        for (j, X) in enumerate([lagrangecxd0(m), duallagrangecxd0(m), raviartthomas(m)])
            areas = [
                CompScienceMeshes.volume(chart(X.geo, i)) for i in 1:numcells(geometry(X))
            ]
            area = sum(areas)

            g, w = H2Trees.adjacencygraph(X)

            b = boundary(geometry(X))

            # edges = BEAST.setminus(skeleton(geometry(X), 1), b))
            edges = BEAST.skeleton(geometry(X), 1)
            faces = BEAST.skeleton(geometry(X), 2)
            Σ = BEAST.connectivity(faces, edges, sign)
            ΣΣ = Σ' * Σ

            println("Testing for mesh ($i / $(length(ms))) for basis function ($j / 3)")
            if isempty(b)
                @test sum(w) ≈ area
            end

            for edge in Graphs.edges(g)
                s, d = src(edge), dst(edge)
                anythingtouches = false
                for sfn in X.fns[s]
                    for dfn in X.fns[d]
                        !iszero(ΣΣ[sfn.cellid, dfn.cellid]) && (anythingtouches = true)
                    end
                end
                @test anythingtouches
            end
        end
    end
end
