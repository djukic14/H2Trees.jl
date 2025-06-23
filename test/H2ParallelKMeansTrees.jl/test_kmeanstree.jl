using CompScienceMeshes, BEAST
using H2Trees
using Test
using SparseArrays
using ParallelKMeans
using LinearAlgebra

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

    tree = H2Trees.KMeansTree(vertices(ms[1]), 10; minvalues=100, n_threads=1)
    tree = H2Trees.KMeansTree(vertices(ms[2]), 10; minvalues=100, n_threads=1)
    tree = H2Trees.KMeansTree(vertices(ms[3]), 10; minvalues=100, n_threads=1)
    tree = H2Trees.KMeansTree(vertices(ms[4]), 10; minvalues=100, n_threads=1)

    points = vertices(ms[1])
    pointsmatrix = reduce(hcat, points)
    numclusters = 10
    kmeansresult = kmeans(pointsmatrix, numclusters; n_threads=1)
    # centers = [SVector{3}(kmeansresult.centers[:, i]) for i in 1:numclusters]
    values = [Vector{Int}() for i in 1:numclusters]
    for i in eachindex(kmeansresult.assignments)
        push!(values[kmeansresult.assignments[i]], i)
    end
    # radii = [0.0 for _ in 1:numclusters]
    # for i in eachindex(radii)
    #     for pointindex in values[i]
    #         dist = norm(centers[i] - points[pointindex])
    #         if dist > radii[i]
    #             radii[i] = dist
    #         end
    #     end
    # end

    # return centers, radii, values
end
