"""
    BoundingBallTree{N,D,T} <: H2ClusterTree

A cluster tree where nodes are bounded by spheres (balls).

This tree structure uses spherical bounding volumes to organize spatial data hierarchically.
Each node is bounded by a ball with a center and radius.

# Type Parameters

  - `N`: The ambient dimension (coordinate space dimension).
  - `D`: The type of nodes.
  - `T`: The type of the radius.

# Fields

  - `nodes::Vector{Node{D}}`: Vector of nodes comprising the tree.
  - `root::Int`: Index of the root node.
  - `center::SVector{N,T}`: Center of the bounding ball of the tree.
  - `radius::T`: Radius of the bounding ball of the tree.
  - `nodesatlevel::Vector{Vector{Int}}`: Vector of vectors, where each inner vector contains the indices of nodes at a specific level.
"""
struct BoundingBallTree{N,D,T} <: H2ClusterTree
    nodes::Vector{Node{D}}
    root::Int
    center::SVector{N,T}
    radius::T
    nodesatlevel::Vector{Vector{Int}}
end

function BoundingBallTree(
    center::C, radius::R; minlevel::Int=1, root::Int=1, balldata=BoundingBallData
) where {N,T<:Number,C<:SVector{N,T},R<:Number}
    rootnode = Node(balldata(Int[], center, radius, minlevel), 0, 0, 0)
    return BoundingBallTree([rootnode], root, center, radius, [Int[]])
end

"""
        BoundingBallTree(points::AbstractVector{SVector{N,T}}, splitwrapper, numsplits::Int; minvalues::Int=numsplits, minlevel::Int=1, root::Int=1, balldata=BoundingBallData, updateradii=boundingsphere, kwargs...)

Construct a `BoundingBallTree` by recursively splitting points into partitions.

This function initializes a root bounding ball over all points and repeatedly calls
`splitwrapper` to partition node values until the stopping criterion is met.

# Arguments

    - `points::AbstractVector{SVector{N,T}}`: Array of points to partition.
    - `splitwrapper`: Callable that returns `(partitions, centers, radii)` for a split.
    - `numsplits::Int`: Number of partitions requested at each split.
    - `minvalues::Int`: Minimum number of points required before splitting a node (default: `numsplits`).
    - `minlevel::Int`: Minimum tree level (default: 1).
    - `root::Int`: Index of root node (default: 1).
    - `balldata`: Data structure for storing bounding ball information (default: `BoundingBallData`).
    - `updateradii`: Function for updating node radii after tree construction (default: `boundingsphere`).
    - `kwargs...`: Additional arguments passed to `splitwrapper`.

# Returns

A `BoundingBallTree` with points organized hierarchically according to `splitwrapper`.

# See also

`KMeansTree`, `MetisTree`, `MetisForest`.
"""
function BoundingBallTree(
    points::AbstractVector{SVector{N,T}},
    splitwrapper,
    numsplits::Int;
    minvalues::Int=numsplits,
    minlevel::Int=1,
    root::Int=1,
    balldata=BoundingBallData,
    updateradii=boundingsphere,
    kwargs...,
) where {N,T}
    center, radius = BoundingSphere.boundingsphere(points)
    tree = BoundingBallTree(center, radius; minlevel=minlevel, root=root, balldata=balldata)
    append!(values(data(tree, root)), collect(1:length(points)))

    _splitboundingballnode!(
        tree,
        points,
        root,
        numsplits,
        splitwrapper;
        minvalues=minvalues,
        balldata=balldata,
        kwargs...,
    )

    _adjustnodesatlevels!(tree)
    updateradii!(tree; update=updateradii)

    return tree
end

function _splitboundingballnode!(
    tree::BoundingBallTree,
    points::AbstractVector{SVector{N,T}},
    node,
    numsplits::Int,
    splitwrapper;
    minvalues::Int=numsplits,
    balldata=BoundingBallData,
    kwargs...,
) where {N,T}
    length(values(tree, node)) <= max(minvalues, numsplits) && return tree

    partitions, centers, radii = splitwrapper(
        points, values(tree, node), numsplits; kwargs...
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
            kwargs...,
        )

        _updatenextsibling!(
            tree, childnode, i == last(eachindex(partitions)) ? 0 : lastnode(tree) + 1
        )
    end
    empty!(values(data(tree, node)))

    return tree
end

function _updatechild!(tree, node, child)
    return tree.nodes[node - root(tree) + 1] = Node(
        data(tree, node), nextsibling(tree, node), parent(tree, node), child
    )
end

function _updatenextsibling!(tree, node, sibling)
    return tree.nodes[node - root(tree) + 1] = Node(
        data(tree, node), sibling, parent(tree, node), firstchild(tree, node)
    )
end

function Base.eltype(::Union{BoundingBallTree{N,D,T},TwoNTree{N,D,T}}) where {N,D,T}
    return SVector{N,T}
end

H2Trees.treetrait(::Type{BoundingBallTree{N,D,T}}) where {N,D,T} = isBoundingBallTree()

# Tis is a very coarse approximation of a bounding sphere.
# See "The Smallest Enclosing Ball of Balls: Combinatorial Structure and Algorithms",
# Fischer (2004) for the right implementation of the SEBB algorithm.
function boundingsphereofspheres(
    center1::A1, radius1::T, center2::A2, radius2::T
) where {T<:Number,A2<:AbstractArray{T},A1<:AbstractArray{T}}
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

Compute a bounding sphere that encloses a tree node and all its children recursively.

This function performs a coarse approximation of the smallest enclosing ball algorithm
by recursively computing bounding spheres for each child node and merging them.

# Arguments

  - `tree`: The bounding ball tree.
  - `node::Int`: Index of the node to compute the bounding sphere for.

# Returns

A tuple `(center, radius)` representing the bounding sphere.
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

Compute a bounding sphere using the maximum child radius (unsafe approximation).

This function computes a bounding sphere by taking the node's center and using the
maximum radius among all child nodes. This is an unsafe approximation that may not
correctly enclose all child nodes and should only be used when the proper bounding
sphere computation has failed or for testing purposes.

!!! warning

    The radii of the tree have not been properly updated. This will lead to incorrect
    results when using the iterators in H2Trees. Proceed at your own risk.

# Arguments

  - `tree`: The bounding ball tree.
  - `node::Int`: Index of the node.

# Returns

A tuple `(center, radius)` where radius is the maximum child radius.
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

Return the node's bounding sphere without updating it.

This function returns the stored center and radius of a node without performing any
computation or update. This can be used when radii updates are intentionally skipped
and existing values should be preserved as-is.

!!! warning

    The radii of the tree have not been updated. This will lead to incorrect results
    when using the iterators in H2Trees. Proceed at your own risk.

# Arguments

  - `tree`: The bounding ball tree.
  - `node::Int`: Index of the node.

# Returns

A tuple `(center, radius)` with the node's current stored values.
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

function updateradii!(tree::BoundingBallTree; update=boundingsphere)
    warning = _updateradiiwarning(update)
    !isnothing(warning) && @warn warning

    for node in DepthFirstIterator(tree)
        center, radius = update(tree, node)
        tree.nodes[node - root(tree) + 1] = Node(
            BoundingBallData(
                values(data(tree, node)),
                SVector(deepcopy(center)),
                radius,
                level(tree, node),
            ),
            nextsibling(tree, node),
            parent(tree, node),
            firstchild(tree, node),
        )
    end
    return tree
end
