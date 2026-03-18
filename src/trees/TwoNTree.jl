"""
    TwoNTree
"""
struct TwoNTree{N,D,T} <: H2ClusterTree
    nodes::Vector{Node{D}}
    root::Int
    center::SVector{N,T}
    halfsize::T
    nodesatlevel::Vector{Vector{Int}}
end

function (tree::H2ClusterTree)(node::Int)
    return tree.nodes[node - H2Trees.root(tree) + 1]
end

function TwoNTree(
    positions,
    minhalfsize;
    minlevel::Int=1,
    root::Int=1,
    minvalues=0,
    maxprotrusion=NaN,
    minsubdividelevel=minlevel + 1,
    computeprotrusion=ComputeProtrusionFunctor(),
)
    rootcenter, rootsize = boundingbox(positions)

    return TwoNTree(
        SVector(rootcenter...),
        positions,
        roothalfsize(rootsize, minhalfsize),
        minhalfsize;
        minlevel=minlevel,
        root=root,
        minvalues=minvalues,
        maxprotrusion=maxprotrusion,
        minsubdividelevel=minsubdividelevel,
        computeprotrusion=computeprotrusion,
    )
end

function TwoNTree(
    center::SVector{N,T}, halfsize::T; minlevel::Int=1, root::Int=1, minvalues=0
) where {N,T}
    rootnode = Node(BoxData(0, Int[], center, halfsize, minlevel), 0, 0, 0)
    return TwoNTree([rootnode], root, center, halfsize, [Int[]])
end

function TwoNTree(
    center::SVector{N,T},
    points::AbstractVector{SVector{N,T}},
    halfsize::T,
    minhalfsize::T;
    minlevel::Int=1,
    root::Int=1,
    minvalues=0,
    maxprotrusion=NaN,
    minsubdividelevel=minlevel + 1,
    computeprotrusion=ComputeProtrusionFunctor(),
) where {N,T}
    tree = TwoNTree(center, halfsize; minlevel=minlevel, root=root)

    addpoints!(
        tree,
        points,
        minhalfsize;
        minsubdividelevel=minsubdividelevel,
        minlevel=minlevel,
        rootsize=halfsize,
        rootcenter=center,
        minvalues=minvalues,
        maxprotrusion=maxprotrusion,
        computeprotrusion=computeprotrusion,
    )
    _adjustnodesatlevels!(tree)

    return tree
end

function sectorcentersize(pt, ct, hs)
    hs = hs / 2
    bl = pt .> ct
    ct = ifelse.(bl, ct .+ hs, ct .- hs)
    sc = sum(b ? 2^(i - 1) : 0 for (i, b) in enumerate(bl))
    return sc, ct, hs
end

# ClusterTrees API #########################################################################

function route!(tree::TwoNTree, state, router)
    point = targetpoint(router)
    smallest_box_size = smallestboxsize(router)

    nodeid, center, size, sfc_state, lvl = state
    size <= smallest_box_size && return state
    !subdivide(router)(router.pointid, lvl) && return state

    target_sector, target_center, target_size = sectorcentersize(point, center, size)
    target_pos = hilbert_positions[sfc_state][target_sector + 1] + 1
    target_sfc_state = hilbert_states[sfc_state][target_sector + 1] + 1
    target_level = lvl + 1

    chds = children(tree, nodeid)
    pos = start(chds)
    while !done(chds, pos)
        child, newpos = next(chds, pos)
        child_sector = sector(data(tree, child))
        child_pos = hilbert_positions[sfc_state][child_sector + 1] + 1
        child_level = level(data(tree, child))
        target_pos < child_pos && break
        if child_sector == target_sector
            return child, target_center, target_size, target_sfc_state, child_level
        end
        pos = newpos
    end

    dat = BoxData(target_sector, Int[], target_center, target_size, target_level)
    child = insert!(chds, dat, pos)

    return child, target_center, target_size, target_sfc_state, target_level
end

function Base.insert!(chd_itr::ChildIterator{<:H2ClusterTree}, item, state)
    prev, next = state
    parent = chd_itr.node

    tree = chd_itr.tree
    push!(tree.nodes, Node(item, next, parent, 0))
    this = lastindex(tree.nodes) + H2Trees.root(tree) - 1
    if prev < 1
        setfirstchild!(tree, parent, this)
    else
        setnextsibling!(tree, prev, this)
    end
    return this
end

# Util functions ###########################################################################

function setfirstchild!(tree::H2ClusterTree, nodeid, child)
    node = tree(nodeid)
    return tree.nodes[nodeid - H2Trees.root(tree) + 1] = Node(
        node.data, node.nextsibling, node.parent, child
    )
end

function setnextsibling!(tree::H2ClusterTree, nodeid, sibling)
    node = tree(nodeid)

    return tree.nodes[nodeid - H2Trees.root(tree) + 1] = Node(
        node.data, sibling, node.parent, node.firstchild
    )
end

"""
    parentcenterminuschildcenter(tree::TwoNTree{N,D,T}, child::Int) where {D,T}

Calculate the difference `r_p-r_c` between the center of the parent `r_p` and the center of
the child node `r_c`.

# Arguments

  - `tree::TwoNTree{N,D,T}`
  - `child::Int`: The index of the child node.

# Returns

  - `SVector{N,T}`: The difference between the center of the parent node and the center of
    the child node.
"""
function parentcenterminuschildcenter(tree::TwoNTree{N,D,T}, child::Int) where {N,D,T}
    return parentcenterminuschildcenter(N, sector(tree, child), halfsize(tree, child))
end

function parentcenterminuschildcenter(N, sector, halfsize)
    ds = digits(sector; base=2, pad=N)
    for i in eachindex(ds)
        ds[i] == 0 && (ds[i] = -1)
    end

    return SVector{N}(ds .* -halfsize)
end

function oppositesector(tree::TwoNTree{N,D,T}, node::Int) where {N,D,T}
    return oppositesector(N, sector(tree, node))
end

function oppositesector(N::Int, sector)
    return (2^N - 1) - sector
end

function comparisonTwoNTree(points, root::Int, roothalfsize, minlevel)
    nlevels = 1
    tree = TwoNTree(points, roothalfsize / 2^nlevels; root=root, minlevel=minlevel)

    while length(tree.nodesatlevel[end]) < length(points)
        nlevels += 1
        tree = TwoNTree(points, roothalfsize / 2^nlevels; root=root, minlevel=minlevel)
    end
    return tree
end

function numberoflevels(halfsize, minhalfsize)
    return ceil(Int, log2(halfsize / minhalfsize))
end

function roothalfsize(rootsize::T, minhalfsize::T) where {T}
    iszero(minhalfsize) && return rootsize
    nlevels = numberoflevels(rootsize, minhalfsize)
    return minhalfsize * T(2)^(nlevels)
end

function addpoint!(
    tree::TwoNTree,
    points,
    pointid,
    smallestboxsize;
    minlevel::Int=H2Trees.level(tree, H2Trees.root(tree)),
    rootsize=H2Trees.halfsize(tree, H2Trees.root(tree)),
    rootcenter=H2Trees.center(tree, H2Trees.root(tree)),
    subdivide=CheckSubdivideFunctor(),
)
    router = Router(smallestboxsize, points, subdivide, pointid)
    root_state = root(tree), rootcenter, rootsize, 1, minlevel
    return update!(tree, root_state, pointid, router) do tree, node, data
        push!(values(tree, node), data)
        node == root(tree) && return nothing
        prnt = parent(tree, node)
        for pointid in values(H2Trees.data(tree, prnt))
            deleteat!(values(H2Trees.data(tree, prnt)), 1)

            addpoint!(
                tree,
                points,
                pointid,
                smallestboxsize;
                minlevel=minlevel,
                rootsize=rootsize,
                rootcenter=rootcenter,
                subdivide=subdivide,
            )
        end
        return nothing
    end
end

function addpoints!(
    tree::TwoNTree,
    points,
    smallestboxsize;
    minlevel::Int=level(tree, root(tree)),
    rootsize=halfsize(tree, root(tree)),
    rootcenter=center(tree, root(tree)),
    minvalues=0,
    minsubdividelevel=minlevel + 1,
    maxprotrusion=NaN,
    computeprotrusion=ComputeProtrusionFunctor(),
)
    subdivide = CheckSubdivideFunctor(
        minvalues,
        maxprotrusion,
        minsubdividelevel,
        computeprotrusion,
        points,
        H2Trees.root(tree),
        minlevel,
        rootsize,
    )
    for i in eachindex(points)
        addpoint!(
            tree,
            points,
            i,
            smallestboxsize;
            minlevel=minlevel,
            rootsize=rootsize,
            rootcenter=rootcenter,
            subdivide=subdivide,
        )
    end
end

"""
    cornerpoints(tree::TwoNTree{N,D,T}, node::Int, i)

Return the corner point of a given node in an N-dimensional TwoNTree.

# Arguments

  - `tree::TwoNTree{N,D,T}`: The tree.
  - `node::Int`: The index of the node.
  - `i`: The corner point index (1 til 2^N).

# Returns

  - The corner point coordinates as a `SVector`.
"""
function cornerpoints(tree::TwoNTree{N,D,T}, node::Int, i) where {N,D,T}
    ds = reverse(digits(i - 1; base=2, pad=N))
    for i in eachindex(ds)
        ds[i] == 0 && (ds[i] = -1)
    end

    return center(tree, node) + SVector{N}(ds .* halfsize(tree, node))
end

H2Trees.treetrait(::Type{TwoNTree{N,D,T}}) where {N,D,T} = isTwoNTree()

"""
    locatepoint(tree, point, level)

Find the node in `tree` at the given `level` whose box contains `point`.

Starting from the root, the function descends the tree by determining which child sector
contains `point` at each level, until the target `level` is reached.

# Arguments

  - `tree`: A `TwoNTree`.
  - `point`: The point to locate.
  - `level`: The tree level at which to find the containing node.

# Returns

  - The node index of the box at `level` that contains `point`.

# Throws

  - An error if no suitable node at `level` is found for `point`.
"""
function locatepoint(tree, point, level)
    tempnode = root(tree)
    for _ in minimumlevel(tree):level
        if H2Trees.level(tree, tempnode) == level
            return tempnode
        end

        newsector, _, _ = sectorcentersize(
            point, center(tree, tempnode), halfsize(tree, tempnode)
        )

        for child in children(tree, tempnode)
            if sector(tree, child) == newsector
                tempnode = child
                break
            end
        end
    end

    return error("No suitable node at level $level found for point $point")
end
