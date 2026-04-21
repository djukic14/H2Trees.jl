function metispartition end # requires Metis.jl to load

"""
        MetisForest(points, graph, weights, numdivisions::Int; splitunconnectedpartitions=false, minlevel::Int=1, minvalues::Int=numdivisions, root::Int=1, balldata=BoundingBallData, updateradii=noboundingsphereupdate)

Construct a `Forest` of `MetisTree`s by partitioning each connected component of a graph.

This function first splits the input graph into connected components, constructs one
`MetisTree` per component, and then remaps local tree indices back to global point indices.

# Arguments

    - `points`: Array of points associated with graph vertices.
    - `graph`: Graph describing connectivity between points.
    - `weights`: Vertex weights used by the METIS partitioning routine.
    - `numdivisions::Int`: Number of partitions requested at each split.
    - `splitunconnectedpartitions`: Whether to further split disconnected METIS partitions (default: `false`).
    - `minlevel::Int`: Minimum tree level (default: 1).
    - `minvalues::Int`: Minimum number of values required before splitting a node (default: `numdivisions`).
    - `root::Int`: Index of root node (default: 1).
    - `balldata`: Data structure for storing bounding ball information (default: `BoundingBallData`).
    - `updateradii`: Function for updating node radii in each component tree (default: `boundingsphere`).

# Returns

A `Forest` containing one `MetisTree` per connected component of `graph`.
"""
function MetisForest(
    points,
    graph,
    weights,
    numdivisions::Int;
    splitunconnectedpartitions=false,
    minlevel::Int=1,
    minvalues::Int=numdivisions,
    root::Int=1,
    balldata=BoundingBallData,
    updateradii=boundingsphere,
)
    trees = []
    for components in connected_components(graph)
        subgraph, localtoglobal = induced_subgraph(graph, components)

        tree = MetisTree(
            points[components],
            subgraph,
            weights[components],
            numdivisions;
            splitunconnectedpartitions=splitunconnectedpartitions,
            minlevel=minlevel,
            minvalues=minvalues,
            root=root,
            balldata=balldata,
            updateradii=updateradii,
        )
        _updatevalues!(tree, localtoglobal)
        push!(trees, tree)
    end
    trees = SVector{length(trees),typeof(trees[begin])}(trees)
    return Forest(trees)
end

"""
    MetisTree(points::AbstractVector{SVector{N,T}}, graph, pointgraphweights, numdivisions::Int; splitunconnectedpartitions=false, minlevel::Int=1, minvalues::Int=numdivisions, root::Int=1, balldata=BoundingBallData, updateradii=noboundingsphereupdate)

Construct a `BoundingBallTree` using METIS graph partitioning to split points.

This function builds a hierarchical tree by recursively partitioning the induced
subgraph of point indices and assigning each partition to a child node.

# Arguments

  - `points::AbstractVector{SVector{N,T}}`: Array of points to partition.
  - `graph`: Graph describing connectivity between points.
  - `pointgraphweights`: Vertex weights used by the METIS partitioning routine.
  - `numdivisions::Int`: Number of partitions requested at each split.
  - `splitunconnectedpartitions`: Whether to further split disconnected METIS partitions (default: `false`).
  - `minvalues::Int`: Minimum number of points required before splitting a node (default: `numdivisions`).
  - `minlevel::Int`: Minimum tree level (default: 1).
  - `root::Int`: Index of root node (default: 1).
  - `balldata`: Data structure for storing bounding ball information (default: `BoundingBallData`).
  - `updateradii`: Function for updating node radii (default: `boundingsphere`).

# Returns

A `BoundingBallTree` with points organized hierarchically using METIS-based splits.
"""
function MetisTree(
    points::AbstractVector{SVector{N,T}},
    graph,
    pointgraphweights,
    numdivisions::Int;
    splitunconnectedpartitions=false,
    minlevel::Int=1,
    minvalues::Int=numdivisions,
    root::Int=1,
    balldata=BoundingBallData,
    updateradii=boundingsphere,
) where {N,T}
    splitwrapper = MetisSplitWrapper(graph, pointgraphweights, splitunconnectedpartitions)
    return BoundingBallTree(
        points,
        splitwrapper,
        numdivisions;
        minvalues=minvalues,
        minlevel=minlevel,
        root=root,
        balldata=balldata,
        updateradii=updateradii,
    )
end

function _updatevalues!(tree, localtoglobal)
    for leaf in leaves(tree)
        globalindices = localtoglobal[values(tree, leaf)]
        empty!(values(tree, leaf))
        append!(values(tree, leaf), globalindices)
    end
    return tree
end

struct MetisSplitWrapper{G,W,S}
    graph::G
    pointgraphweights::W
    splitunconnectedpartitions::S
end

function (f::MetisSplitWrapper)(points, globalpointids, numsplits; kwargs...)
    return metiswrapper(
        f.graph,
        points,
        globalpointids,
        f.pointgraphweights,
        numsplits;
        splitunconnectedpartitions=f.splitunconnectedpartitions,
        kwargs...,
    )
end

function metiswrapper(
    graph::Graphs.SimpleGraph,
    points::Vector{SVector{N,T}},
    globalpointids::Vector{Int},
    pointgraphweights,
    numdivisions::Int;
    kwargs...,
) where {N,T}
    subgraph, _ = induced_subgraph(graph, globalpointids)
    part = metispartition(
        subgraph, pointgraphweights[globalpointids], numdivisions; kwargs...
    )
    partitions = [Vector{Int}() for _ in 1:maximum(part)]
    for (i, p) in enumerate(part)
        push!(partitions[p], globalpointids[i])
    end

    filter!(!isempty, partitions)
    centers = Vector{SVector{N,T}}(undef, length(partitions))
    radii = Vector{T}(undef, length(partitions))

    for i in eachindex(partitions)
        centers[i], radii[i] = BoundingSphere.boundingsphere(points[partitions[i]])
    end

    return partitions, centers, radii
end
