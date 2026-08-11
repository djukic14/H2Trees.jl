"""
    DepthFirstIterator(tree, node=root(tree))

Iterate over `node` and all of its descendants in depth-first order.
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

"""
    _LeafFunctor(tree)

Callable predicate returning `isleaf(tree, node)`.
"""
struct _LeafFunctor{T}
    tree::T
end

function (f::_LeafFunctor)(node::Int)
    return isleaf(f.tree, node)
end

"""
    leaves(tree, node::Int)

Return the leaf nodes below `node`.

When `node` is the root, the cached [`TreeIndex`](@ref) leaf list is reused and
copied for caller ownership.
"""
function leaves(tree, node::Int=H2Trees.root(tree))
    node == H2Trees.root(tree) && return copy(treeindex(tree).leaves)
    return collect(
        Int, Iterators.filter(_LeafFunctor(tree), DepthFirstIterator(tree, node))
    )
end

"""
    ChildIterator(tree, node)

Iterate over the direct children of `node`.
"""
struct ChildIterator{T,N<:Integer}
    tree::T
    node::N
end

Base.IteratorSize(cv::ChildIterator) = Base.SizeUnknown()
Base.eltype(::ChildIterator{T,N}) where {T,N} = N

"""
    ParentUpwardsIterator(tree, node)

Iterate from `parent(tree, node)` upward until the root has been reached.
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
    return prnt, prnt
end

function Base.iterate(itr::ParentUpwardsIterator, node)
    prnt = parent(itr.tree, node)
    if iszero(prnt)
        return nothing
    end

    return prnt, prnt
end

# Utils DepthFirstIterator #################################################################

"""
    NodeInformation(info)

Typed wrapper for a child-iterator state, including `nothing`.
"""
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

"""
    StackElement(childreniterator, information)

Stack frame used by [`DepthFirstIterator`](@ref).
"""
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

# Buckets `ids` by `sectorcentersize`'s sector (always `0:2^N-1`, see `_uniformseparationdepth`
# for why a plain `Vector` + `isassigned` beats a `Dict` here). Callers reconstruct a bucket's
# child center on demand via `sectorcenter(sector, center, childhalfsize)` instead of this
# function also allocating/returning a `centers` vector alongside the buckets; `sectorcenter` is
# a handful of arithmetic ops, cheaper than the allocation it replaces.
"""
    _bucketbysector(points, ids, center, halfsize, nsectors)

Group point ids by their child sector relative to `center` and `halfsize`.

Only assigned buckets are initialized. Callers reconstruct child centers from
the sector id when needed.
"""
function _bucketbysector(points, ids, center, halfsize, nsectors)
    buckets = Vector{Vector{Int}}(undef, nsectors)
    # A crude but cheap even split: better than starting every bucket from empty and repeatedly
    # doubling as it grows, without tracking exact per-sector counts up front.
    buckethint = max(1, length(ids) ÷ nsectors)
    for i in ids
        sc, _, _ = sectorcentersize(points[i], center, halfsize)
        idx = sc + 1
        if !isassigned(buckets, idx)
            buckets[idx] = sizehint!(Int[], buckethint)
        end
        push!(buckets[idx], i)
    end
    return buckets
end

"""
    _anyprotrudes(points, ids, center, halfsize, maxprotrusion, computeprotrusion)

Return whether any point id protrudes beyond the accepted protrusion threshold.
"""
function _anyprotrudes(points, ids, center, halfsize, maxprotrusion, computeprotrusion)
    for id in ids
        computeprotrusion(center, halfsize, id) >= maxprotrusion && return true
    end
    return false
end

"""
    _validateprotrusion(protrusion)

Validate the protrusion-check object used by bulk `TwoNTree` construction.
"""
function _validateprotrusion(protrusion)
    protrusion isa NoProtrusionCheck ||
        protrusion isa ProtrusionCheck ||
        throw(
            ArgumentError(
                "protrusion must be a `NoProtrusionCheck` or `ProtrusionCheck`, got " *
                "$(typeof(protrusion))",
            ),
        )
    return nothing
end
