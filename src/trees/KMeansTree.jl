function kmeanswrapper end # requires ParallelKMeans.jl to load

struct KMeansSplitWrapper end

function (f::KMeansSplitWrapper)(points, globalpointids, numsplits; kwargs...)
    return kmeanswrapper(points, globalpointids, numsplits; kwargs...)
end

"""
    KMeansTree(points::AbstractVector{SVector{N,T}}, numberofclusters::Int; minvalues::Int=numberofclusters, minlevel::Int=1, root::Int=1, balldata=BoundingBallData, updateradii=boundingsphere, kwargs...)

Construct a `BoundingBallTree` using k-means clustering to partition the given points.

This function recursively clusters the points using k-means algorithm to build a hierarchical
tree structure. Each node in the tree contains a cluster of points bounded by a sphere.

# Arguments

  - `points::AbstractVector{SVector{N,T}}`: Array of points to partition.
  - `numberofclusters::Int`: Number of clusters for k-means at each split.
  - `minvalues::Int`: Minimum number of points required before splitting a node (default: `numberofclusters`).
  - `minlevel::Int`: Minimum tree level (default: 1).
  - `root::Int`: Index of root node (default: 1).
  - `balldata`: Data structure for storing bounding ball information (default: `BoundingBallData`).
  - `updateradii`: Function for updating node radii (default: `boundingsphere`).
  - `kwargs...`: Additional arguments passed to the k-means wrapper function.

# Returns

A `BoundingBallTree` with points organized hierarchically by k-means clustering.
"""
function KMeansTree(
    points::AbstractVector{SVector{N,T}},
    numberofclusters::Int;
    minvalues::Int=numberofclusters,
    minlevel::Int=1,
    root::Int=1,
    balldata=BoundingBallData,
    updateradii=boundingsphere,
    kwargs...,
) where {N,T}
    splitwrapper = KMeansSplitWrapper()
    return BoundingBallTree(
        points,
        splitwrapper,
        numberofclusters;
        minvalues=minvalues,
        minlevel=minlevel,
        root=root,
        balldata=balldata,
        updateradii=updateradii,
        kwargs...,
    )
end
