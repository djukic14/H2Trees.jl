module H2ParallelKMeansTrees
using ParallelKMeans
using StaticArrays, LinearAlgebra
using H2Trees
import H2Trees: kmeanswrapper

function _computeradius(center, points)
    maxdist = zero(eltype(center))
    for point in points
        dist = norm(center - point)
        if dist > maxdist
            maxdist = dist
        end
    end
    return maxdist
end

function kmeanswrapper(
    points::Vector{SVector{N,T}},
    globalpointids::Vector{Int},
    numberofclusters::Int;
    kwargs...,
) where {N,T}
    pointsmatrix = reduce(hcat, points[globalpointids])
    kresult = kmeans(pointsmatrix, numberofclusters; kwargs...)
    centers = [SVector{N}(kresult.centers[:, i]) for i in axes(kresult.centers, 2)]

    partitions = [Vector{Int}() for _ in 1:maximum(kresult.assignments)]
    for (i, p) in enumerate(kresult.assignments)
        push!(partitions[p], globalpointids[i])
    end

    radii = [_computeradius(centers[i], points[partitions[i]]) for i in eachindex(centers)]

    return partitions, centers, radii
end

end # module H2ParallelKMeansTrees
