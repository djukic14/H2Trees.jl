using CompScienceMeshes, BEAST
using H2Trees
using Test
using SparseArrays
using ParallelKMeans
using LinearAlgebra
using PlotlyJS
@testset "KMeansTree" begin
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

    for i in eachindex(ms)
        tree = H2Trees.KMeansTree(vertices(ms[i]), 10; minvalues=100, n_threads=1)
        @test H2Trees.testwellseparatedness(tree)
    end
end
