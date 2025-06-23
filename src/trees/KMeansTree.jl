function KMeansTree(
    points::AbstractVector{SVector{N,T}},
    numberofclusters::Int;
    minvalues::Int=numberofclusters,
    minlevel::Int=1,
    root::Int=1,
    balldata=BoundingBallData,
    kwargs...,
) where {N,T}
    pointsmatrix = reduce(hcat, points)
    kmeansresult = kmeans(pointsmatrix, 1; kwargs...)
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
    kresult = kmeans(pointsmatrix, numberofclusters; kwargs...)
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
