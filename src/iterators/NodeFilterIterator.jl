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
    NodeFilterIterator(tree, node::Int, filter)

Iterate over same-level nodes that pass `filter`.

For ordinary trees, candidates come from `tree` at `level(tree, node)` and the
predicate is called as `filter(tree, candidate, node)`. For block trees, `node`
is interpreted as a trial-tree node and candidates come from the test tree.
"""
function NodeFilterIterator(tree, node::Int, filter)
    return NodeFilterIterator(tree, node::Int, treetrait(tree), filter)
end

function NodeFilterIterator(tree, node::Int, ::AbstractTreeTrait, filter)
    return Iterators.filter(
        _LevelFilterFunctor(tree, node, filter), LevelIterator(tree, level(tree, node))
    )
end

function NodeFilterIterator(tree, trialnode::Int, ::isBlockTree, filter)
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
