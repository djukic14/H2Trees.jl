"""
    BoundingBallTree{N,D,T} <: H2ClusterTree

Cluster tree whose nodes are bounded by balls.

`N` is the ambient dimension, `D` is the node-data type, and `T` is the radius
type. The cached [`TreeIndex`](@ref) is stored in a `Ref` so
`rebuildtreeindex!` can replace it after topology-changing operations.
"""
struct BoundingBallTree{N,D,T,I} <: H2ClusterTree
    nodes::Vector{Node{D}}
    root::Int
    center::SVector{N,T}
    radius::T
    index::I
end

"""
    BoundingBallTree(nodes, root, center, radius, nodesatlevel)

Construct a `BoundingBallTree` from already materialized nodes and level
storage, then rebuild the cached tree index.
"""
function BoundingBallTree(
    nodes::Vector{Node{D}},
    root::Int,
    center::SVector{N,T},
    radius::T,
    nodesatlevel::Vector{Vector{Int}},
) where {N,D,T}
    minlevel = isempty(nodesatlevel) ? 0 : level(first(nodes).data)
    maxlevel = isempty(nodesatlevel) ? 0 : minlevel + length(nodesatlevel) - 1
    index = Ref(TreeIndex(nodesatlevel, Int[], Int[], minlevel, maxlevel))
    tree = BoundingBallTree{N,D,T,Base.RefValue{TreeIndex}}(
        nodes, root, center, radius, index
    )
    return rebuildtreeindex!(tree)
end

function BoundingBallTree{N,D,T}(
    nodes::Vector{Node{D}},
    root::Int,
    center::SVector{N,T},
    radius::T,
    nodesatlevel::Vector{Vector{Int}},
) where {N,D,T}
    return BoundingBallTree(nodes, root, center, radius, nodesatlevel)
end

"""
    BoundingBallTree(center, radius; minlevel=1, root=1, balldata=BoundingBallData)

Construct a one-node bounding-ball tree.
"""
function BoundingBallTree(
    center::C, radius::R; minlevel::Int=1, root::Int=1, balldata=BoundingBallData
) where {N,T<:Number,C<:SVector{N,T},R<:Number}
    rootnode = Node(balldata(Int[], center, radius, minlevel), 0, 0, 0)
    return BoundingBallTree(
        [rootnode],
        root,
        center,
        radius,
        Ref(TreeIndex([[root]], [root], [root], minlevel, minlevel)),
    )
end

"""
    BoundingBallTree(points::AbstractVector{SVector{N,T}}; builder::BoundingBallTreeBuilder)

Construct a `BoundingBallTree` from points.

Points are recursively split with `builder.splitter` until the builder's
stopping criteria are met. Prefer the canonical [`buildtree`](@ref) entry point;
this method is what it forwards to.
"""
function BoundingBallTree(
    points::AbstractVector{SVector{N,T}}; builder::BoundingBallTreeBuilder
) where {N,T}
    center, radius = BoundingSphere.boundingsphere(points)
    minlevel = _resolve_minlevel(builder.minlevel, 1)
    tree = BoundingBallTree(
        center, radius; minlevel=minlevel, root=builder.root, balldata=builder.balldata
    )
    append!(values(data(tree, builder.root)), 1:length(points))

    _splitboundingballnode!(
        tree,
        points,
        builder.root,
        builder.numsplits,
        builder.splitter;
        minvalues=builder.minvalues,
        balldata=builder.balldata,
        maxlevel=builder.maxlevel,
        builder.splitterkwargs...,
    )

    rebuildtreeindex!(tree)
    updateradii!(tree; update=builder.updateradii)

    return tree
end

"""
    _splitboundingballnode!(tree, points, node, numsplits, splitwrapper; kwargs...)

Recursively split one bounding-ball node and append its children to `tree`.

The split wrapper returns child value partitions, centers, and radii. Internal
node values are cleared after their children are created, so stored values live
on leaves.
"""
function _splitboundingballnode!(
    tree::BoundingBallTree,
    points::AbstractVector{SVector{N,T}},
    node,
    numsplits::Int,
    splitwrapper;
    minvalues::Int=numsplits,
    balldata=BoundingBallData,
    maxlevel=typemax(Int),
    kwargs...,
) where {N,T}
    nodevalues = values(data(tree, node))
    length(nodevalues) <= max(minvalues, numsplits) && return tree
    level(tree, node) >= maxlevel && return tree
    partitions, centers, radii = _callsplitwrapper(
        splitwrapper, points, nodevalues, level(tree, node), numsplits; kwargs...
    )
    _updatechild!(tree, node, lastnode(tree) + 1)

    for i in eachindex(partitions)
        dat = balldata(partitions[i], centers[i], radii[i], level(tree, node) + 1)
        childnode = lastnode(tree) + 1

        push!(tree.nodes, Node(dat, 0, node, 0))
        _splitboundingballnode!(
            tree,
            points,
            childnode,
            numsplits,
            splitwrapper;
            minvalues=minvalues,
            balldata=balldata,
            maxlevel=maxlevel,
            kwargs...,
        )

        _updatenextsibling!(
            tree, childnode, i == last(eachindex(partitions)) ? 0 : lastnode(tree) + 1
        )
    end

    empty!(values(data(tree, node)))

    return tree
end

"""
    _callsplitwrapper(splitwrapper, points, globalpointids, level, numsplits; kwargs...)

Call a bounding-ball split wrapper.

Split wrappers may accept either `(points, values, level, numsplits; kwargs...)`
or the older `(points, values, numsplits; kwargs...)` shape.
"""
function _callsplitwrapper(
    splitwrapper, points, globalpointids, level, numsplits; kwargs...
)
    if applicable(splitwrapper, points, globalpointids, level, numsplits)
        return splitwrapper(points, globalpointids, level, numsplits; kwargs...)
    elseif applicable(splitwrapper, points, globalpointids, numsplits)
        return splitwrapper(points, globalpointids, numsplits; kwargs...)
    end

    return throw(
        ArgumentError(
            "splitwrapper must support either (points, values, level, numsplits; kwargs...) " *
            "or (points, values, numsplits; kwargs...)",
        ),
    )
end

"""
    _copyatlevel(data::BoundingBallData, level)

Copy node data while replacing its stored level.
"""
function _copyatlevel(data::BoundingBallData, level::Int)
    return typeof(data)(copy(values(data)), center(data), radius(data), level)
end

"""
    balanceleaves!(tree::BoundingBallTree)

Extend shallower leaves with unary child nodes until all leaves have the same level.

The added child nodes keep the same bounding ball and values as the original leaf. The
original leaf becomes an internal node and its values are cleared, matching the usual tree
layout where values live on leaves.
"""
function balanceleaves!(tree::BoundingBallTree)
    leafnodes = leaves(tree)
    isempty(leafnodes) && return tree

    leaflevel = maximum(level.(Ref(tree), leafnodes))
    for leaf in leafnodes
        currentnode = leaf
        while level(tree, currentnode) < leaflevel
            childnode = lastnode(tree) + 1
            childdata = _copyatlevel(data(tree, currentnode), level(tree, currentnode) + 1)
            push!(tree.nodes, Node(childdata, 0, currentnode, 0))
            _updatechild!(tree, currentnode, childnode)
            empty!(values(data(tree, currentnode)))
            currentnode = childnode
        end
    end

    rebuildtreeindex!(tree)
    return tree
end

"""
    _updatechild!(tree, node, child)

Set `child` as the first child of `node`.
"""
function _updatechild!(tree, node, child)
    return tree.nodes[node - root(tree) + 1] = Node(
        data(tree, node), nextsibling(tree, node), parent(tree, node), child
    )
end

"""
    _updatenextsibling!(tree, node, sibling)

Set `sibling` as the next sibling of `node`.
"""
function _updatenextsibling!(tree, node, sibling)
    return tree.nodes[node - root(tree) + 1] = Node(
        data(tree, node), sibling, parent(tree, node), firstchild(tree, node)
    )
end

"""
    eltype(tree::Union{BoundingBallTree,TwoNTree})

Return the static-vector coordinate type stored by tree geometry.
"""
function Base.eltype(::Union{BoundingBallTree{N,D,T},TwoNTree{N,D,T}}) where {N,D,T}
    return SVector{N,T}
end

"""
    treetrait(::Type{<:BoundingBallTree})

Return [`isBoundingBallTree`](@ref).
"""
H2Trees.treetrait(::Type{<:BoundingBallTree}) = isBoundingBallTree()

"""
    boundingsphereofspheres(center1, radius1, center2, radius2)

Return a ball enclosing two input balls.

This is a coarse two-ball merge, not a full smallest-enclosing-ball-of-balls
algorithm.
"""
function boundingsphereofspheres(
    center1::A1, radius1::T, center2::A2, radius2::T
) where {T<:Number,A2<:AbstractArray{T},A1<:AbstractArray{T}}
    # This is a very coarse approximation of a bounding sphere.
    # See "The Smallest Enclosing Ball of Balls: Combinatorial Structure and Algorithms",
    # Fischer (2004) for the right implementation of the SEBB algorithm.
    difference = center1 - center2
    differencenorm = norm(difference)

    if (differencenorm + radius2 <= radius1)
        # ball2 is inside ball1
        return center1, radius1

    elseif (differencenorm + radius1 <= radius2)
        # ball1 is inside ball2
        return center2, radius2
    else
        center =
            T(0.5) * (center1 + center2 + (radius1 - radius2) * difference / differencenorm)
        radius = T(0.5) * (radius1 + radius2 + differencenorm)

        return center, radius
    end
end

"""
    boundingsphere(tree, node::Int)

Compute a bounding sphere for `node` from its children.

The implementation recursively merges child spheres with
[`boundingsphereofspheres`](@ref), so it is a conservative approximation rather
than an exact smallest enclosing ball.
"""
function boundingsphere(tree, node::Int)
    centerbuffer = similar(center(tree, node))
    centerbuffer .= center(tree, node)
    rds = radius(tree, node)
    for (i, child) in enumerate(children(tree, node))
        i == 1 && (centerbuffer .= center(tree, child))

        centerbuffer, rds = boundingsphereofspheres(
            centerbuffer, rds, center(tree, child), radius(tree, child)
        )
    end

    return centerbuffer, rds
end

"""
    unsafemaxradiusboundingsphere(tree, node::Int)

Return `center(tree, node)` and the maximum radius among `node` and its children.

!!! warning

    This can fail to enclose all child balls. It is intended only as an unsafe
    fallback or testing hook.
"""
function unsafemaxradiusboundingsphere(tree, node::Int)
    rds = radius(tree, node)
    for child in children(tree, node)
        rds = max(rds, radius(tree, child))
    end

    return center(tree, node), rds
end

"""
    noboundingsphereupdate(tree, node::Int)

Return the stored center and radius of `node`.

!!! warning

    Use only when stale radii are intentional; near/far iterators rely on
    correct radii.
"""
function noboundingsphereupdate(tree, node::Int)
    return center(tree, node), radius(tree, node)
end

function _updateradiiwarning(::typeof(unsafemaxradiusboundingsphere))
    return "Radii of the tree have not been properly updated. This will lead to incorrect results when using the iterators in H2Trees.  Proceed at your own risk."
end
function _updateradiiwarning(::typeof(noboundingsphereupdate))
    return "Radii of the tree have not been updated. This will lead to incorrect results when using the iterators in H2Trees.  Proceed at your own risk."
end
_updateradiiwarning(_) = nothing

"""
    updateradii!(tree::BoundingBallTree; update=boundingsphere)

Update every node radius and center with `update(tree, node)`.

`update` may be [`boundingsphere`](@ref), [`noboundingsphereupdate`](@ref), or a
custom function returning `(center, radius)`.
"""
function updateradii!(tree::BoundingBallTree; update=boundingsphere)
    warning = _updateradiiwarning(update)
    !isnothing(warning) && @warn warning

    for node in DepthFirstIterator(tree)
        center, radius = update(tree, node)
        tree.nodes[node - root(tree) + 1] = Node(
            BoundingBallData(
                values(data(tree, node)), SVector(center), radius, level(tree, node)
            ),
            nextsibling(tree, node),
            parent(tree, node),
            firstchild(tree, node),
        )
    end
    return tree
end
