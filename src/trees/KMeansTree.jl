"""
    kmeanswrapper(points, globalpointids, numberofclusters; kwargs...)

Extension hook implemented by `H2ParallelKMeansTrees`.

The core package declares the function so [`KMeansTree`](@ref) can be defined
without loading `ParallelKMeans.jl`.
"""
function kmeanswrapper end

"""
    KMeansSplitWrapper()

Split-wrapper adapter used by [`KMeansTree`](@ref)'s underlying
[`BoundingBallTree`](@ref) construction.
"""
struct KMeansSplitWrapper end

"""
    (wrapper::KMeansSplitWrapper)(points, globalpointids, level, numsplits; kwargs...)

Call [`kmeanswrapper`](@ref) for the selected point ids.

`level` is accepted to match the generic bounding-ball split-wrapper protocol;
the K-means splitter itself does not use it.
"""
function (f::KMeansSplitWrapper)(points, globalpointids, level, numsplits; kwargs...)
    return kmeanswrapper(points, globalpointids, numsplits; kwargs...)
end

"""
    KMeansTree(points; builder=KMeansTreeBuilder())

Build a ball tree by recursively splitting point sets with K-means.

`KMeansTree` is implemented as a [`BoundingBallTree`](@ref) with a
K-means-backed split wrapper. Load the `H2ParallelKMeansTrees` extension, or a
package that activates it, so [`kmeanswrapper`](@ref) is available.
"""
function KMeansTree(
    points::AbstractVector{SVector{N,T}}; builder::KMeansTreeBuilder=KMeansTreeBuilder()
) where {N,T}
    splitwrapper = KMeansSplitWrapper()
    ballbuilder = BoundingBallTreeBuilder(;
        splitter=splitwrapper,
        numsplits=builder.numberofclusters,
        minvalues=builder.minvalues,
        minlevel=builder.minlevel,
        root=builder.root,
        balldata=builder.balldata,
        updateradii=builder.updateradii,
        splitterkwargs=builder.splitterkwargs,
    )
    return BoundingBallTree(points; builder=ballbuilder)
end
