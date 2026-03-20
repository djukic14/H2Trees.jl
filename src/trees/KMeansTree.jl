function kmeanswrapper end # requires ParallelKMeans.jl to load

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
    pointsmatrix = reduce(hcat, points)
    kmeansresult = kmeanswrapper(pointsmatrix, 1; kwargs...)
    center = SVector{N}(kmeansresult.centers[:, 1])
    radius = _computeradius(center, points)

    tree = BoundingBallTree(center, radius; minlevel=minlevel, root=root, balldata=balldata)
    append!(values(data(tree, root)), collect(1:length(points)))

    splitnode!(
        tree,
        points,
        root,
        numberofclusters;
        minvalues=minvalues,
        balldata=balldata,
        kwargs...,
    )

    _adjustnodesatlevels!(tree)
    updateradii!(tree; update=updateradii)

    return tree
end

function splitnode!(
    tree::BoundingBallTree,
    points::AbstractVector{SVector{N,T}},
    node,
    numberofclusters::Int;
    minvalues::Int=numberofclusters,
    balldata=BoundingBallData,
    kwargs...,
) where {N,T}
    length(values(tree, node)) < max(minvalues, numberofclusters) && return tree

    pointsmatrix = reduce(hcat, points[values(tree, node)])
    kresult = kmeanswrapper(pointsmatrix, numberofclusters; kwargs...)
    centers = [SVector{N}(kresult.centers[:, i]) for i in axes(kresult.centers, 2)]
    vals = [Vector{Int}() for _ in 1:numberofclusters]
    for i in eachindex(kresult.assignments)
        push!(vals[kresult.assignments[i]], values(tree, node)[i])
    end
    radii = [_computeradius(centers[i], points[vals[i]]) for i in eachindex(centers)]

    _updatechild!(tree, node, lastnode(tree) + 1)

    for i in eachindex(centers)
        dat = balldata(vals[i], centers[i], radii[i], level(tree, node) + 1)
        childnode = lastnode(tree) + 1
        push!(tree.nodes, Node(dat, 0, node, 0))
        splitnode!(
            tree,
            points,
            childnode,
            numberofclusters;
            minvalues=minvalues,
            balldata=balldata,
            kwargs...,
        )

        _updatenextsibling!(
            tree, childnode, i == last(eachindex(centers)) ? 0 : lastnode(tree) + 1
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

# Tis is a very coarse approximation of a bounding sphere.
# See "The Smallest Enclosing Ball of Balls: Combinatorial Structure and Algorithms",
# Fischer (2004) for the right implementation of the SEBB algorithm.
function boundingsphere(
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

        centerbuffer, rds = boundingsphere(
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
    @warn "Radii of the tree have not been properly updated. This will lead to incorrect results when using the iterators in H2Trees.  Proceed at your own risk."
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
    @warn "Radii of the tree have not been updated. This will lead to incorrect results when using the iterators in H2Trees.  Proceed at your own risk."
    return center(tree, node), radius(tree, node)
end

function updateradii!(tree::BoundingBallTree; update=boundingsphere)
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
