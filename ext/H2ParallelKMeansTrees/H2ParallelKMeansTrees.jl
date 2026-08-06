module H2ParallelKMeansTrees
using ParallelKMeans
using StaticArrays, LinearAlgebra
using H2Trees
import H2Trees: kmeanswrapper

"""
    _computeradius(center, points, ids)

Return the maximum distance from `center` to the points indexed by `ids`.
"""
function _computeradius(center, points, ids)
    maxdist = zero(eltype(center))
    for i in ids
        dist = norm(center - points[i])
        if dist > maxdist
            maxdist = dist
        end
    end
    return maxdist
end

"""
    kmeanswrapper(points, globalpointids, numberofclusters; kwargs...)

Partition a subset of points with `ParallelKMeans.kmeans`.

`globalpointids` selects the points to cluster and is preserved in the returned
partitions. The result is `(partitions, centers, radii)`, with empty clusters
removed before returning.
"""
function kmeanswrapper(
    points::Vector{SVector{N,T}},
    globalpointids::Vector{Int},
    numberofclusters::Int;
    kwargs...,
) where {N,T}
    pointsmatrix = Matrix{T}(undef, N, length(globalpointids))
    for (j, i) in enumerate(globalpointids)
        @inbounds for d in 1:N
            pointsmatrix[d, j] = points[i][d]
        end
    end
    kresult = kmeans(pointsmatrix, numberofclusters; kwargs...)
    # K-means can return empty requested clusters; allocate for every returned
    # center first, then filter empty partitions before constructing tree nodes.
    nclusters = size(kresult.centers, 2)
    allcenters = [SVector{N}(view(kresult.centers, :, i)) for i in 1:nclusters]

    rawpartitions = [Vector{Int}() for _ in 1:nclusters]
    for (i, p) in enumerate(kresult.assignments)
        push!(rawpartitions[p], globalpointids[i])
    end

    # Keep `partitions` and `centers` single-assignment so the `radii`
    # comprehension stays type-stable instead of boxing captured variables.
    nonempty = findall(!isempty, rawpartitions)
    partitions = length(nonempty) == nclusters ? rawpartitions : rawpartitions[nonempty]
    centers = length(nonempty) == nclusters ? allcenters : allcenters[nonempty]
    radii = [_computeradius(centers[i], points, partitions[i]) for i in eachindex(centers)]

    return partitions, centers, radii
end

end # module H2ParallelKMeansTrees
