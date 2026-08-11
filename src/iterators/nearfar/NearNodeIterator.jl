
"""
    NearNodeIterator(tree, node::Int; isnear=isnear, nearlists=nothing)

Iterate over same-level nodes in `tree` that are near `node`.

The predicate is called as `isnear(tree, candidate, node)`.

`nearlists` optionally supplies a [`NearListCache`](@ref) to draw candidates from
instead of scanning the level; see [`nearlistcache`](@ref) for when that is valid.
"""
function NearNodeIterator(tree, node::Int; isnear=isnear, nearlists=nothing)
    return NodeFilterIterator(tree, node, isnear; nearlists=nearlists)
end

"""
    NearNodeIterator(testtree, trialtree, trialnode::Int; isnear=isnear)

Iterate over nodes in `testtree` near `trialnode` in `trialtree`.

Only nodes at `level(trialtree, trialnode)` are considered. The predicate is
called as `isnear(testtree, trialtree, candidate, trialnode)`.
"""
function NearNodeIterator(testtree, trialtree, trialnode::Int; isnear=isnear)
    return NearNodesAtAnchorLevel(testtree, trialtree, trialnode; isnear=isnear)
end

function NearNodesAtAnchorLevel(iteratedtree, anchortree, anchornode::Int; isnear=isnear)
    return SameLevelFilteredIterator(iteratedtree, anchortree, anchornode, isnear)
end

"""
    FarNodeIterator(tree, node::Int; isfar=isfar)

Iterate over same-level nodes in `tree` that are far from `node`.
"""
function FarNodeIterator(tree, node::Int; isfar=isfar)
    return NodeFilterIterator(tree, node, isfar)
end

"""
    FarNodeIterator(testtree, trialtree, trialnode::Int; isfar=isfar)

Iterate over same-level nodes in `testtree` that are far from `trialnode` in
`trialtree`.
"""
function FarNodeIterator(testtree, trialtree, trialnode::Int; isfar=isfar)
    return FarNodesAtAnchorLevel(testtree, trialtree, trialnode; isfar=isfar)
end

function FarNodesAtAnchorLevel(iteratedtree, anchortree, anchornode::Int; isfar=isfar)
    return SameLevelFilteredIterator(iteratedtree, anchortree, anchornode, isfar)
end

struct _LeafPredicateFunctor{P}
    predicate::P
end

function (f::_LeafPredicateFunctor)(tree, nodea, nodeb)
    return isleaf(tree, nodea) && f.predicate(tree, nodea, nodeb)
end

function (f::_LeafPredicateFunctor)(testtree, trialtree, testnode, trialnode)
    return isleaf(testtree, testnode) &&
           f.predicate(testtree, trialtree, testnode, trialnode)
end

function _getindicesstorage(::Val{:flattened})
    return Int[]
end

function _storevalues!(indices, tree, node, ::Val{:flattened})
    return appendvalues!(indices, tree, node)
end

function _appendnodevalues!(
    indices, iteratedtree, iterator, predicate, storevalues, anchorargs...
)
    for node in iterator(anchorargs...; predicate=predicate)
        _storevalues!(indices, iteratedtree, node, storevalues)
    end

    anchortree = anchorargs[end - 1]
    anchornode = anchorargs[end]

    leafpredicate = _LeafPredicateFunctor(predicate)
    for parent in ParentUpwardsIterator(anchortree, anchornode)
        for node in iterator(anchorargs[1:(end - 1)]..., parent; predicate=leafpredicate)
            _storevalues!(indices, iteratedtree, node, storevalues)
        end
    end
    return indices
end

# Carries the near cache without teaching the shared append helper about near-only state.
struct _NearNodeIteratorFactory{L}
    nearlists::L
end

function (f::_NearNodeIteratorFactory)(tree, node; predicate)
    return NearNodeIterator(tree, node; isnear=predicate, nearlists=f.nearlists)
end

function (f::_NearNodeIteratorFactory)(testtree, trialtree, trialnode; predicate)
    return NearNodeIterator(testtree, trialtree, trialnode; isnear=predicate)
end

function _farnodeiterator(tree, node; predicate)
    return FarNodeIterator(tree, node; isfar=predicate)
end

function _farnodeiterator(testtree, trialtree, trialnode; predicate)
    return FarNodeIterator(testtree, trialtree, trialnode; isfar=predicate)
end

"""
    nearnodevalues(tree, node::Int; isnear=isnear, storevalues=Val{:flattened}(), nearlists=nothing)

Collect values stored in near nodes for `node` in a single tree.

The traversal first visits near nodes at the level of `node`, then walks upward
through parents and includes near leaf nodes. By default, returns a flattened
`Vector{Int}`.

Pass `nearlists` when calling this for many nodes of the same tree; see
[`nearlistcache`](@ref).
"""
function nearnodevalues(
    tree, node::Int; isnear=isnear, storevalues=Val{:flattened}(), nearlists=nothing
)
    indices = _getindicesstorage(storevalues)
    appendnearnodevalues!(
        indices, tree, node; isnear=isnear, storevalues=storevalues, nearlists=nearlists
    )
    return indices
end

"""
    appendnearnodevalues!(indices, tree, node::Int; isnear=isnear, storevalues=Val{:flattened}(), nearlists=nothing)

Append single-tree near-node values to `indices` and return `indices`.
"""
function appendnearnodevalues!(
    indices,
    tree,
    node::Int;
    isnear=isnear,
    storevalues=Val{:flattened}(),
    nearlists=nothing,
)
    return _appendnodevalues!(
        indices, tree, _NearNodeIteratorFactory(nearlists), isnear, storevalues, tree, node
    )
end

"""
    nearnodevalues(testtree, trialtree, trialnode::Int; isnear=isnear, storevalues=Val{:flattened}())

Collect values from near nodes in `testtree` for a reference node in `trialtree`.

The traversal first visits near nodes relative to `trialnode`, then walks upward
through trial-tree parents and adds near leaf nodes. By default, returns a
flattened `Vector{Int}`.
"""
function nearnodevalues(
    testtree, trialtree, trialnode::Int; isnear=isnear, storevalues=Val{:flattened}()
)
    indices = _getindicesstorage(storevalues)
    appendnearnodevalues!(
        indices, testtree, trialtree, trialnode; isnear=isnear, storevalues=storevalues
    )
    return indices
end

"""
    appendnearnodevalues!(indices, testtree, trialtree, trialnode::Int; isnear=isnear, storevalues=Val{:flattened}())

Append two-tree near-node values from `testtree` to `indices` and return
`indices`.
"""
function appendnearnodevalues!(
    indices,
    testtree,
    trialtree,
    trialnode::Int;
    isnear=isnear,
    storevalues=Val{:flattened}(),
)
    return _appendnodevalues!(
        indices,
        testtree,
        _NearNodeIteratorFactory(nothing),
        isnear,
        storevalues,
        testtree,
        trialtree,
        trialnode,
    )
end

"""
    farnodevalues(tree, node::Int; isfar=isfar)

Collect values stored in far nodes for `node` in a single tree.
"""
function farnodevalues(tree, node::Int; isfar=isfar)
    indices = Int[]
    appendfarnodevalues!(indices, tree, node; isfar=isfar)
    return indices
end

"""
    appendfarnodevalues!(indices, tree, node::Int; isfar=isfar)

Append single-tree far-node values to `indices` and return `indices`.
"""
function appendfarnodevalues!(indices, tree, node::Int; isfar=isfar)
    return _appendnodevalues!(
        indices, tree, _farnodeiterator, isfar, Val(:flattened), tree, node
    )
end

"""
    farnodevalues(testtree, trialtree, trialnode::Int; isfar=isfar)

Collect values from far nodes in `testtree` for a reference node in `trialtree`.
"""
function farnodevalues(testtree, trialtree, trialnode::Int; isfar=isfar)
    indices = Int[]
    appendfarnodevalues!(indices, testtree, trialtree, trialnode; isfar=isfar)
    return indices
end

"""
    appendfarnodevalues!(indices, testtree, trialtree, trialnode::Int; isfar=isfar)

Append two-tree far-node values from `testtree` to `indices` and return
`indices`.
"""
function appendfarnodevalues!(indices, testtree, trialtree, trialnode::Int; isfar=isfar)
    return _appendnodevalues!(
        indices,
        testtree,
        _farnodeiterator,
        isfar,
        Val(:flattened),
        testtree,
        trialtree,
        trialnode,
    )
end
