struct _LevelFilterFunctor{F,T,TN}
    tree::T
    trialnode::TN
    filter::F
end

function (f::_LevelFilterFunctor)(x::Int)
    return f.filter(f.tree, x, f.trialnode)
end

struct _LevelBlockTreeFilterFunctor{F,TE,TR,TN}
    testtree::TE
    trialtree::TR
    trialnode::TN
    filter::F
end

function (f::_LevelBlockTreeFilterFunctor)(testnode::Int)
    return f.filter(f.testtree, f.trialtree, testnode, f.trialnode)
end

"""
    NodeFilterIterator(tree, node::Int, filter)

Returns an iterator that returns nodes at the same level as `node`, for which the function
`filter(tree, nodea, nodeb)` or the function `filter(testtree, trialtree, testnode, trialnode)`
return true. In the case of a tree with the trait `isBlockTree` it is assumed that node is
belonging to the `trialtree`.

# Arguments

  - `tree`: The tree to iterate over
  - `node::Int`: The node to start from
  - `filter`: A function that takes a tree and two nodes and returns a boolean

# Returns

An iterator over nodes at the same level as `node` that pass the filter
"""
function NodeFilterIterator(tree, node::Int, filter)
    return NodeFilterIterator(tree, node::Int, treetrait(tree), filter)
end

# for any tree that is not a BlockTree
function NodeFilterIterator(tree, node::Int, ::AbstractTreeTrait, filter)
    return Iterators.filter(
        _LevelFilterFunctor(tree, node, filter), LevelIterator(tree, level(tree, node))
    )
end

# a BlockTree has a testtree and a trialtree
# we assume that the node passed in this case is from the trialtree
# the nodes returned are from the testtree
function NodeFilterIterator(tree, trialnode::Int, ::isBlockTree, filter)
    return NodeFilterIterator(testtree(tree), trialtree(tree), trialnode, filter)
end

"""
    NodeFilterIterator(testtree, trialtree, trialnode::Int, filter)

Returns an iterator over nodes at the same level as `trialnode` in `trialtree`,
for which the function `filter(testtree, trialtree, testnode, trialnode)` return true.

In the case of a tree with the trait `isBlockTree`, it is assumed that `trialnode` is
belonging to the `trialtree`.

# Arguments

  - `testtree`: The tree to iterate over
  - `trialtree`: The tree to which the `trialnode` belongs to
  - `trialnode::Int`: The node to start from
  - `filter`: A function that takes a test tree, a trial tree, a test node, and a trial node and returns a boolean

# Returns

An iterator over nodes at the same level as `trialnode` that pass the filter
"""
function NodeFilterIterator(testtree, trialtree, trialnode::Int, filter)
    return Iterators.filter(
        _LevelBlockTreeFilterFunctor(testtree, trialtree, trialnode, filter),
        LevelIterator(testtree, level(trialtree, trialnode)),
    )
end
