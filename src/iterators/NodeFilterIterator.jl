struct LevelFilterFunctor{F,T,TN}
    tree::T
    trialnode::TN
    filter::F
end

function (f::LevelFilterFunctor)(x::Int)
    return f.filter(f.tree, x, f.trialnode)
end

struct LevelBlockTreeFilterFunctor{F,TE,TR,TN}
    testtree::TE
    trialtree::TR
    trialnode::TN
    filter::F
end

function (f::LevelBlockTreeFilterFunctor)(testnode::Int)
    return f.filter(f.testtree, f.trialtree, testnode, f.trialnode)
end

function NodeFilterIterator(tree, trialnode::Int, filter)
    return NodeFilterIterator(tree, trialnode::Int, treetrait(tree), filter)
end

# for any tree that is not a BlockTree
function NodeFilterIterator(tree, trialnode::Int, ::H2Trees.AbstractTreeTrait, filter)
    nodelevel = H2Trees.level(tree, trialnode)
    return Iterators.filter(
        LevelFilterFunctor(tree, trialnode, filter), LevelIterator(tree, nodelevel)
    )
end

# a BlockTree has a testtree and a trialtree
# we assume that the node passed in this case is from the trialtree
# the nodes returned are from the testtree
function NodeFilterIterator(tree, trialnode::Int, ::isBlockTree, filter)
    return NodeFilterIterator(testtree(tree), trialtree(tree), trialnode, filter)
end

function NodeFilterIterator(testtree, trialtree, trialnode::Int, filter)
    triallevel = level(trialtree, trialnode)
    return Iterators.filter(
        LevelBlockTreeFilterFunctor(testtree, trialtree, trialnode, filter),
        LevelIterator(testtree, triallevel),
    )
end
