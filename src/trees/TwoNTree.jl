struct TwoNTree{N,D,T,P} <: H2ClusterTree
    nodes::Vector{Node{D}}
    root::Int
    center::SVector{N,T}
    halfsize::T
    nodesatlevel::Vector{Vector{Int}}
    pointsunique::P
end

function (tree::H2ClusterTree)(node::Int)
    return tree.nodes[node - H2Trees.root(tree) + 1]
end

function boxdata(sector, values, center, halfsize, level)
    return BoxData(sector, values, center, halfsize, level)
end

function TwoNTree(
    positions,
    minhalfsize;
    minlevel::Int=1,
    root::Int=1,
    boxdata=boxdata,
    pointsunique=H2Trees.UniquePoints(),
    minvalues=0,
)
    rootcenter, rootsize = boundingbox(positions)

    nlevels = numberoflevels(rootsize, minhalfsize)

    return TwoNTree(
        SVector(rootcenter...),
        positions,
        roothalfsize(minhalfsize, nlevels),
        minhalfsize;
        minlevel=minlevel,
        root=root,
        boxdata=boxdata,
        pointsunique=pointsunique,
        minvalues=minvalues,
    )
end

function TwoNTree(
    center::SVector{N,T},
    halfsize::T;
    minlevel::Int=1,
    root::Int=1,
    boxdata=boxdata,
    pointsunique=H2Trees.UniquePoints(),
    minvalues=0,
) where {N,T}
    rootnode = Node(boxdata(0, Int[], center, halfsize, minlevel), 0, 0, 0)
    return TwoNTree([rootnode], root, center, halfsize, [Int[]], pointsunique)
end

function TwoNTree(
    center::SVector{N,T},
    points::AbstractArray{SVector{N,T},1},
    halfsize::T,
    minhalfsize::T;
    minlevel::Int=1,
    root::Int=1,
    boxdata=boxdata,
    pointsunique=H2Trees.UniquePoints(),
    minvalues=0,
) where {N,T}
    tree = TwoNTree(
        center,
        halfsize;
        minlevel=minlevel,
        root=root,
        boxdata=boxdata,
        pointsunique=pointsunique,
    )

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

function route!(tree::TwoNTree, state, destination; boxdata=boxdata)
    point = destination.target_point
    smallest_box_size = destination.smallest_box_size
    minvalues = destination.minvalues

    nodeid, center, size, sfc_state, lvl = state
    size <= smallest_box_size && return state
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
        node.data, node.next_sibling, node.parent, child
    )
end

function setnextsibling!(tree::H2ClusterTree, nodeid, sibling)
    node = tree(nodeid)

    return tree.nodes[nodeid - H2Trees.root(tree) + 1] = Node(
        node.data, sibling, node.parent, node.first_child
    )
end

"""
    parentcenterminuschildcenter(tree::TwoNTree{3,D,T}, child::Int) where {D,T}

Calculate the difference `r_p-r_c` between the center of the parent `r_p` and the center of
the child node `r_c`. This currently only works for octrees.

# Arguments

  - `tree::TwoNTree{3,D,T}`
  - `child::Int`: The index of the child node.

# Returns

  - `SVector{3,T}`: The difference between the center of the parent node and the center of
    the child node.
"""
function parentcenterminuschildcenter(tree::TwoNTree{3,D,T}, child::Int) where {D,T}
    return parentcenterminuschildcenter(sector(tree, child), halfsize(tree, child))
end

"""
    parentcenterminuschildcenter(sector::Int, halfsizechild::T) where {T}

Calculate the difference `r_p-r_c` between the center of the parent `r_p` and the center of
the child node `r_c`. This currently only works for octrees.

# Arguments

  - `sector::Int`: The sector of the child.
  - `halfsizechild::T`: Half size of the child box.

# Returns

  - `SVector{3,T}`: The difference between the parent center and the child center.

# Example

```
julia> parentcenterminuschildcenter(0, 0.5)
3-element SVector{3, Float64}:
 0.5
 0.5
 0.5
```
"""
function parentcenterminuschildcenter(sector::Int, halfsizechild::T) where {T}
    @match sector begin
        0 => return SVector{3,T}(halfsizechild, halfsizechild, halfsizechild)
        1 => return SVector{3,T}(-halfsizechild, halfsizechild, halfsizechild)
        2 => return SVector{3,T}(halfsizechild, -halfsizechild, halfsizechild)
        3 => return SVector{3,T}(-halfsizechild, -halfsizechild, halfsizechild)
        4 => return SVector{3,T}(halfsizechild, halfsizechild, -halfsizechild)
        5 => return SVector{3,T}(-halfsizechild, halfsizechild, -halfsizechild)
        6 => return SVector{3,T}(halfsizechild, -halfsizechild, -halfsizechild)
        7 => return SVector{3,T}(-halfsizechild, -halfsizechild, -halfsizechild)
    end
end

function oppositesector(tree::TwoNTree, node::Int)
    return oppositesector(sector(tree, node))
end

function oppositesector(sector::Int)
    @match sector begin
        0 => return 7
        1 => return 6
        2 => return 5
        3 => return 4
        4 => return 3
        5 => return 2
        6 => return 1
        7 => return 0
    end
end

function numberoflevels(halfsize, minhalfsize)
    return ceil(Int, log2(halfsize / minhalfsize))
end

function roothalfsize(minhalfsize::T, nlevels) where {T}
    return minhalfsize * T(2)^(nlevels)
end

function addpoint!(
    tree::TwoNTree,
    point,
    pointid,
    smallestboxsize;
    minlevel::Int=H2Trees.level(tree, H2Trees.root(tree)),
    rootsize=H2Trees.halfsize(tree, H2Trees.root(tree)),
    rootcenter=H2Trees.center(tree, H2Trees.root(tree)),
    boxdata=boxdata,
    minvalues=0,
)
    router = Router(smallestboxsize, point, minvalues)
    root_state = root(tree), rootcenter, rootsize, 1, minlevel
    return update!(tree, root_state, pointid, router; boxdata=boxdata) do tree, node, data
        return push!(values(tree, node), data)
    end
end

function addpoints!(
    tree::TwoNTree,
    points,
    smallestboxsize;
    minlevel::Int=level(tree, root(tree)),
    rootsize=halfsize(tree, root(tree)),
    rootcenter=center(tree, root(tree)),
    boxdata=boxdata,
    minvalues=0,
)
    for i in eachindex(points)
        addpoint!(
            tree,
            points[i],
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
    cornerpoints(tree::TwoNTree{2,D,T}, node::Int, i)

Return the corner point of a given node in a two-dimensional quadtree.

# Arguments

  - `tree::TwoNTree{2,D,T}`: The quadtree.
  - `node::Int`: The index of the node.
  - `i`: The corner point index (1, 2, 3, or 4).

# Returns

  - The corner point coordinates as a `SVector`.
"""
function cornerpoints(tree::TwoNTree{2,D,T}, node::Int, i) where {D,T}
    center = H2Trees.center(tree, node)
    halfsize = H2Trees.halfsize(tree, node)
    @assert i in 1:4

    @match i begin
        1 => return center + SVector(-halfsize, -halfsize)
        2 => return center + SVector(-halfsize, halfsize)
        3 => return center + SVector(halfsize, -halfsize)
        4 => return center + SVector(halfsize, halfsize)
    end
end

"""
    cornerpoints(tree::TwoNTree{3,D,T}, node::Int, i)

Return the corner point of a given node in a three-dimensional octree.

# Arguments

  - `tree::TwoNTree{3,D,T}`: The quadtree.
  - `node::Int`: The index of the node.
  - `i`: The corner point index (1 til 8).

# Returns

  - The corner point coordinates as a `SVector`.
"""
function cornerpoints(tree::TwoNTree{3,D,T}, node::Int, i) where {D,T}
    center = H2Trees.center(tree, node)
    halfsize = H2Trees.halfsize(tree, node)
    @assert i in 1:8

    @match i begin
        1 => return center + SVector(-halfsize, -halfsize, -halfsize)
        2 => return center + SVector(-halfsize, -halfsize, halfsize)
        3 => return center + SVector(-halfsize, halfsize, -halfsize)
        4 => return center + SVector(-halfsize, halfsize, halfsize)
        5 => return center + SVector(halfsize, -halfsize, -halfsize)
        6 => return center + SVector(halfsize, -halfsize, halfsize)
        7 => return center + SVector(halfsize, halfsize, -halfsize)
        8 => return center + SVector(halfsize, halfsize, halfsize)
    end
    #TODO: instead of writing one function for 2 and 3 dimensions --> write a function that converts i into bits and use the bit representation to calculate the corner point
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

H2Trees.treetrait(::Type{TwoNTree{N,D,T,P}}) where {N,D,T,P} = isTwoNTree()

function uniquepointstree(tree::TwoNTree{N,D,T,H2Trees.NonUniquePoints}) where {N,D,T}
    return false
end
