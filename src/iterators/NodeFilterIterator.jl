struct _LevelFilterFunctor{F,T,TN}
    tree::T
    trialnode::TN
    filter::F
end

function (f::_LevelFilterFunctor)(x::Int)
    return f.filter(f.tree, x, f.trialnode)
end

struct _LevelBlockTreeFilterFunctor{F,TE,TR,TN}
    iteratedtree::TE
    anchortree::TR
    anchornode::TN
    filter::F
end

function (f::_LevelBlockTreeFilterFunctor)(iteratednode::Int)
    return f.filter(f.iteratedtree, f.anchortree, iteratednode, f.anchornode)
end

"""
    NodeFilterIterator(tree, node::Int, filter; nearlists=nothing)

Iterate over same-level nodes that pass `filter`.

For ordinary trees, candidates come from `tree` at `level(tree, node)` and the
predicate is called as `filter(tree, candidate, node)`. For block trees, `node`
is interpreted as a trial-tree node and candidates come from the test tree.

Passing a [`NearListCache`](@ref) as `nearlists` narrows the candidate set to that
node's cached near list instead of scanning the whole level. That is only correct
when the cache contains every node `filter` could accept, such as near-type filters
built by [`nearlistcache`](@ref). It is never correct for a far filter.
"""
function NodeFilterIterator(tree, node::Int, filter; nearlists=nothing)
    return NodeFilterIterator(tree, node::Int, treetrait(tree), filter, nearlists)
end

function NodeFilterIterator(tree, node::Int, ::AbstractTreeTrait, filter, nearlists=nothing)
    return Iterators.filter(
        _LevelFilterFunctor(tree, node, filter), _nearcandidates(nearlists, tree, node)
    )
end

# No `nearlists` here: `nearlistcache` declines block trees, so it is always `nothing`.
function NodeFilterIterator(tree, trialnode::Int, ::isBlockTree, filter, nearlists=nothing)
    return SameLevelFilteredIterator(testtree(tree), trialtree(tree), trialnode, filter)
end

"""
    NodeFilterIterator(testtree, trialtree, trialnode::Int, filter)

Iterate over nodes in `testtree` at `level(trialtree, trialnode)`.

The predicate is called as `filter(testtree, trialtree, candidate, trialnode)`.
"""
function NodeFilterIterator(testtree, trialtree, trialnode::Int, filter)
    return SameLevelFilteredIterator(testtree, trialtree, trialnode, filter)
end

"""
    SameLevelFilteredIterator(iteratedtree, anchortree, anchornode::Int, predicate)

Return nodes from `iteratedtree` at the same level as `anchornode` in `anchortree`
for which `predicate(iteratedtree, anchortree, iteratednode, anchornode)` is true.
"""
function SameLevelFilteredIterator(iteratedtree, anchortree, anchornode::Int, predicate)
    return Iterators.filter(
        _LevelBlockTreeFilterFunctor(iteratedtree, anchortree, anchornode, predicate),
        LevelIterator(iteratedtree, level(anchortree, anchornode)),
    )
end
