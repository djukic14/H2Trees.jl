import Base.insert!

function update!(f, tree, state, data, target; kwargs...)
    while true
        next_state = route!(tree, state, target; kwargs...)
        next_state == state && break
        state = next_state
    end
    node = first(state)
    f(tree, node, data)
    return node
end

"""
    DepthFirstIterator{T,N<:Integer}

Traverses the tree in a depth first manner. If no node is specified the tree is
traversed from the root node.

# Fields

  - `tree::T`: The tree.
  - `node::Int`: The node from which the tree is traversed.

# Methods

## `DepthFirstIterator(tree, node)`

Creates a new `DepthFirstIterator` instance, traversing the tree from the specified `node`.

## `DepthFirstIterator(tree)`

Creates a new `DepthFirstIterator` instance, traversing the tree from the root node.
"""
struct DepthFirstIterator{T,N<:Integer}
    tree::T
    node::N
end

function DepthFirstIterator(tree)
    return DepthFirstIterator(tree, root(tree))
end

Base.IteratorSize(::DepthFirstIterator) = Base.SizeUnknown()
Base.eltype(::DepthFirstIterator{T,N}) where {T,N} = N

function Base.iterate(itr::DepthFirstIterator)
    chitr = children(itr.tree, itr.node)
    stack = [StackElement(chitr, iterate(chitr))]
    return iterate(itr, stack)
end

function Base.iterate(itr::DepthFirstIterator, stack)
    isempty(stack) && return nothing
    while true
        info = information(last(stack))
        if !isnothing(info)
            n = node(info)
            chitr = children(itr.tree, n)
            push!(stack, eltype(stack)(chitr, iterate(chitr)))
        else
            pop!(stack)
            isempty(stack) && return (itr.node, stack)
            chitr = childreniterator(last(stack))
            info = information(last(stack))
            n, s = node(info), state(info)
            stack[end] = eltype(stack)(chitr, iterate(chitr, s))
            return n, stack
        end
    end
end

struct LeafFunctor{T}
    tree::T
end

function (f::LeafFunctor)(node::Int)
    return H2Trees.isleaf(f.tree, node)
end

"""
    leaves(tree, node::Int)

Returns an iterator over the leaf nodes in the tree, starting from the specified `node`. If
no node is specified the tree is traversed from the root node.

# Arguments

  - `tree`: The tree to search for leaf nodes.
  - `node::Int`: The node from which to start the search.

# Returns

An iterator over the leaf nodes in the tree.
"""
function leaves(tree, node::Int=H2Trees.root(tree))
    return collect(
        Int, Iterators.filter(LeafFunctor(tree), H2Trees.DepthFirstIterator(tree, node))
    )
end

"""
    ChildIterator{T,N<:Integer}

An iterator over the children of a node in a tree.

# Fields

  - `tree::T`: The tree.
  - `node::N`: The node whose children are being iterated over.
"""
struct ChildIterator{T,N<:Integer}
    tree::T
    node::N
end

Base.IteratorSize(cv::ChildIterator) = Base.SizeUnknown()
Base.eltype(::ChildIterator{T,N}) where {T,N} = N

"""
    ParentUpwardsIterator{T,N<:Int}

ParentUpwardsIterator is an iterator that iterates over all parent nodes of a given node in
a tree until the root is reached. The last node is the node 0.

# Fields

  - `tree::T`: The tree.
  - `node::Int`: The node over which parents is iterated.
"""
struct ParentUpwardsIterator{T,N<:Integer}
    tree::T
    node::N
end

Base.IteratorSize(::ParentUpwardsIterator) = Base.SizeUnknown()
Base.eltype(::ParentUpwardsIterator{T,N}) where {T,N} = N

function Base.iterate(itr::ParentUpwardsIterator{T,N}) where {T,N}
    if itr.node == root(itr.tree)
        return nothing
    end

    prnt = parent(itr.tree, itr.node)

    stack = N[prnt]
    return prnt, stack
end

function Base.iterate(itr::ParentUpwardsIterator, stack)
    isempty(stack) && return nothing

    node = stack[begin]
    popfirst!(stack)
    prnt = parent(itr.tree, node)
    if iszero(prnt)
        return nothing
    else
        pushfirst!(stack, prnt)
    end

    return prnt, stack
end

# Utils DepthFirstIterator #################################################################

struct NodeInformation{N}
    info::Union{Nothing,N}
    function NodeInformation(info)
        return new{typeof(info)}(info)
    end
    function NodeInformation{N}(::Nothing) where {N}
        return new{N}(nothing)
    end
    function NodeInformation{N}(info::N) where {N}
        return new{N}(info)
    end
end

node(next::NodeInformation) = next.info[1]

state(next::NodeInformation) = next.info[2]

Base.isnothing(x::NodeInformation) = isnothing(x.info)

struct StackElement{C,N}
    chitr::C
    info::NodeInformation{N}

    function StackElement(chitr, info)
        return new{typeof(chitr),typeof(info)}(chitr, NodeInformation(info))
    end

    function StackElement{C,N}(chitr, info) where {C,N}
        return new{C,N}(chitr, NodeInformation{N}(info))
    end
end

childreniterator(s::StackElement) = s.chitr

information(s::StackElement) = s.info

#TODO: hilbert_positions and hilbert_states for N≠3
struct Router{T,P,M}
    smallest_box_size::T
    target_points::P
    pointid::Int
    minvalues::M
end

function targetpoint(router::Router)
    return router.target_points[router.pointid]
end

function minvalues(router::Router)
    return router.minvalues
end

function smallestboxsize(router::Router)
    return router.smallest_box_size
end

const hilbert_states = [
    [1, 2, 3, 2, 4, 5, 3, 5],
    [2, 6, 0, 7, 8, 8, 0, 7],
    [0, 9, 10, 9, 1, 1, 11, 11],
    [6, 0, 6, 11, 9, 0, 9, 8],
    [11, 11, 0, 7, 5, 9, 0, 7],
    [4, 4, 8, 8, 0, 6, 10, 6],
    [5, 7, 5, 3, 1, 1, 11, 11],
    [6, 1, 6, 10, 9, 4, 9, 10],
    [10, 3, 1, 1, 10, 3, 5, 9],
    [4, 4, 8, 8, 2, 7, 2, 3],
    [7, 2, 11, 2, 7, 5, 8, 5],
    [10, 3, 2, 6, 10, 3, 4, 4],
]

const hilbert_positions = [
    [0, 1, 3, 2, 7, 6, 4, 5],
    [0, 7, 1, 6, 3, 4, 2, 5],
    [0, 3, 7, 4, 1, 2, 6, 5],
    [2, 3, 1, 0, 5, 4, 6, 7],
    [4, 3, 5, 2, 7, 0, 6, 1],
    [6, 5, 1, 2, 7, 4, 0, 3],
    [4, 7, 3, 0, 5, 6, 2, 1],
    [6, 7, 5, 4, 1, 0, 2, 3],
    [2, 5, 3, 4, 1, 6, 0, 7],
    [2, 1, 5, 6, 3, 0, 4, 7],
    [4, 5, 7, 6, 3, 2, 0, 1],
    [6, 1, 7, 0, 5, 2, 4, 3],
]

function start(itr::ChildIterator{<:H2ClusterTree})
    return (0, firstchild(itr.tree, itr.node))
end

function done(itr::ChildIterator{<:H2ClusterTree}, state)
    _, this = state
    this < 1 && return true
    sibling_par = parent(itr.tree, this)
    sibling_par != itr.node && return true
    return false
end

function next(itr::ChildIterator{<:H2ClusterTree}, state)
    prev, this = state
    nxt = nextsibling(itr.tree, this)
    return (this, (this, nxt))
end

function Base.iterate(itr::ChildIterator{<:H2ClusterTree}, st=start(itr))
    return done(itr, st) ? nothing : next(itr, st)
end
