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
    positions, minhalfsize; minlevel::Int=1, root::Int=1, boxdata=BoxData, minvalues=0
)
    rootcenter, rootsize = boundingbox(positions)

    return TwoNTree(
        SVector(rootcenter...),
        positions,
        roothalfsize(rootsize, minhalfsize),
        minhalfsize;
        minlevel=minlevel,
        root=root,
        boxdata=boxdata,
        minvalues=minvalues,
    )
end

function TwoNTree(
    center::SVector{N,T},
    halfsize::T;
    minlevel::Int=1,
    root::Int=1,
    boxdata=BoxData,
    minvalues=0,
) where {N,T}
    rootnode = Node(boxdata(0, Int[], center, halfsize, minlevel), 0, 0, 0)
    return TwoNTree([rootnode], root, center, halfsize, [Int[]])
end

function TwoNTree(
    center::SVector{N,T},
    points::AbstractVector{SVector{N,T}},
    halfsize::T,
    minhalfsize::T;
    minlevel::Int=1,
    root::Int=1,
    boxdata=BoxData,
    minvalues=0,
) where {N,T}
    tree = TwoNTree(center, halfsize; minlevel=minlevel, root=root, boxdata=boxdata)

    addpoints!(
        tree,
        points,
        minhalfsize;
        minlevel=minlevel,
        rootsize=halfsize,
        rootcenter=center,
        boxdata=boxdata,
        minvalues=minvalues,
    )
    _adjustnodesatlevels!(tree)

    return tree
end

# ClusterTrees API #########################################################################

function route!(tree::TwoNTree, state, router; boxdata=BoxData)
    point = targetpoint(router)
    smallest_box_size = smallestboxsize(router)
    minvals = minvalues(router)

    nodeid, center, size, sfc_state, lvl = state
    size <= smallest_box_size && return state
    isleaf(tree, nodeid) && length(values(tree, nodeid)) < minvals && return state

    target_sector, target_center, target_size = sector_center_size(point, center, size)
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

    dat = boxdata(target_sector, Int[], target_center, target_size, target_level)
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
    boxdata=BoxData,
    minvalues=0,
)
    router = Router(smallestboxsize, points, pointid, minvalues)
    root_state = root(tree), rootcenter, rootsize, 1, minlevel
    return update!(tree, root_state, pointid, router; boxdata=boxdata) do tree, node, data
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
                boxdata=boxdata,
                minvalues=minvalues,
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
    boxdata=BoxData,
    minvalues=0,
)
    for i in eachindex(points)
        addpoint!(
            tree,
            points,
            i,
            smallestboxsize,
            ;
            minlevel=minlevel,
            rootsize=rootsize,
            rootcenter=rootcenter,
            boxdata=boxdata,
            minvalues=minvalues,
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

# More API for H2Trees.jl ##################################################################

struct LeafFunctor{T}
    tree::T
end

function (f::LeafFunctor)(node::Int)
    return H2Trees.isleaf(f.tree, node)
end

function leaves(tree::H2ClusterTree, node::Int)
    return collect(
        Int, Iterators.filter(LeafFunctor(tree), H2Trees.DepthFirstIterator(tree, node))
    )
end

H2Trees.treetrait(::Type{TwoNTree{N,D,T}}) where {N,D,T} = isTwoNTree()
