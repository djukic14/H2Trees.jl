
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
        )
        _updatevalues!(tree, localtoglobal)
        push!(trees, tree)
    end
    trees = SVector{length(trees),typeof(trees[begin])}(trees)
    return Forest(trees)
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

# uses Metis for split
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
    updateradii=noboundingsphereupdate,
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
    centers = Vector{SVector{N,T}}(undef, length(partitions))
    radii = Vector{T}(undef, length(partitions))

    for i in eachindex(partitions)
        centers[i], radii[i] = BoundingSphere.boundingsphere(points[partitions[i]])
    end

    return partitions, centers, radii
end
