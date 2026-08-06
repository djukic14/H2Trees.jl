"""
    metispartition(graph, vertexweights, numberofdivisions; kwargs...)

Extension hook implemented by `H2MetisTrees`.

The core package declares the function so [`MetisTree`](@ref) and
[`MetisForest`](@ref) can be defined without loading `Metis.jl`.
"""
function metispartition end

"""
    MetisForest(points, graph, weights; builder::MetisForestBuilder)

Construct a `Forest` of `MetisTree`s from the connected components of `graph`.

Each component is built as its own [`MetisTree`](@ref), then leaf values are
remapped from component-local indices back to the original global point ids.
Prefer the canonical [`buildforest`](@ref) entry point.
"""
function MetisForest(
    points, graph, weights; builder::MetisForestBuilder=MetisForestBuilder()
)
    trees = []
    for components in connected_components(graph)
        subgraph, localtoglobal = induced_subgraph(graph, components)

        tree = MetisTree(
            points[components], subgraph, weights[components]; builder=builder.treebuilder
        )
        _updatevalues!(tree, localtoglobal)
        push!(trees, tree)
    end
    trees = SVector{length(trees),typeof(trees[begin])}(trees)
    return Forest(trees)
end

"""
    MetisTree(points::AbstractVector{SVector{N,T}}, graph, pointgraphweights; builder::MetisTreeBuilder)

Construct a ball tree using METIS graph partitioning.

The tree is implemented as a [`BoundingBallTree`](@ref) whose split wrapper
recursively partitions the induced subgraph of each node's point ids. Prefer the
canonical [`buildtree`](@ref) entry point.
"""
function MetisTree(
    points::AbstractVector{SVector{N,T}},
    graph,
    pointgraphweights;
    builder::MetisTreeBuilder=MetisTreeBuilder(),
) where {N,T}
    splitwrapper = MetisSplitWrapper(
        graph, pointgraphweights, builder.splitunconnectedpartitions
    )
    ballbuilder = BoundingBallTreeBuilder(;
        splitter=splitwrapper,
        numsplits=builder.numdivisions,
        minvalues=builder.minvalues,
        minlevel=builder.minlevel,
        maxlevel=builder.maxlevel,
        root=builder.root,
        balldata=builder.balldata,
        updateradii=builder.updateradii,
    )
    return BoundingBallTree(points; builder=ballbuilder)
end

"""
    _updatevalues!(tree, localtoglobal)

Remap leaf values from component-local indices to global point ids.

Used after constructing one component of a [`MetisForest`](@ref).
"""
function _updatevalues!(tree, localtoglobal)
    for leaf in leaves(tree)
        leafvalues = values(tree, leaf)
        for i in eachindex(leafvalues)
            leafvalues[i] = localtoglobal[leafvalues[i]]
        end
    end
    return tree
end

"""
    MetisSplitWrapper(graph, pointgraphweights, splitunconnectedpartitions)

Split-wrapper state used by [`MetisTree`](@ref)'s underlying
[`BoundingBallTree`](@ref) construction.
"""
struct MetisSplitWrapper{G,W,S}
    graph::G
    pointgraphweights::W
    splitunconnectedpartitions::S
end

"""
    (wrapper::MetisSplitWrapper)(points, globalpointids, level, numsplits; kwargs...)

Partition the subgraph induced by `globalpointids`.

`level` is accepted to match the generic bounding-ball split-wrapper protocol;
the METIS splitter itself does not use it.
"""
function (f::MetisSplitWrapper)(points, globalpointids, level, numsplits; kwargs...)
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

"""
    metiswrapper(graph, points, globalpointids, pointgraphweights, numdivisions; kwargs...)

Build child partitions, centers, and radii for one METIS tree node.

The returned point ids remain global to the original point vector. Empty
partitions are removed before centers and bounding radii are computed.
"""
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
