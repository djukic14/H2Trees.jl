struct _WellSeparatedIteratorFunctor{IN}
    iswellseparated::IN
end

function (f::_WellSeparatedIteratorFunctor)(tree)
    return WellSeparatedIterator(
        tree, treetrait(tree); iswellseparated=f.iswellseparated(tree)
    )
end

"""
    WellSeparatedIterator(; isnear=nothing, iswellseparated=nothing)

Constructs a functor that returns a `WellSeparatedIterator` if provided a tree.
Two nodes are considered well-separated if their parents are near each other and
the nodes themselves are far apart.
This assumes that child clusters are completely inside their parent clusters.

# Arguments

  - `isnear`: a function that takes a tree as input and returns another function. This returned function is then used to evaluate the `isnear` criterion.
  - `iswellseparated`: a function that takes a tree as input and returns another function. This returned function is then used to evaluate the `iswellseparated` criterion.

# Returns

A `WellSeparatedIteratorFunctor` that returns a `WellSeparatedIterator` if provided with a tree.

# Throws

  - `error`: if both `isnear` and `iswellseparated` are provided, or if neither is provided.
"""
function WellSeparatedIterator(;
    iswellseparated=nothing, isnear=isnothing(iswellseparated) ? H2Trees.isnear() : nothing
)
    if !((isnear !== nothing) ⊻ (iswellseparated !== nothing))
        error("Supply one of (not both) isnear or iswellseparated")
    end

    filter = if isnothing(iswellseparated)
        H2Trees.iswellseparated(; isnear=isnear)
    else
        iswellseparated
    end

    return _WellSeparatedIteratorFunctor(filter)
end

struct _WellSeparatedIteratorNotBlockTreeFunctor{IN}
    iswellseparated::IN
end

function (f::_WellSeparatedIteratorNotBlockTreeFunctor)(tree, node)
    return WellSeparatedIterator(tree, node; iswellseparated=f.iswellseparated)
end

function WellSeparatedIterator(tree, ::AbstractTreeTrait; iswellseparated=iswellseparated)
    return _WellSeparatedIteratorNotBlockTreeFunctor(iswellseparated)
end

struct _WellSeparatedIteratorBlockTreeFunctor{IN}
    iswellseparated::IN
end

function (f::_WellSeparatedIteratorBlockTreeFunctor)(testtree, trialtree, trialnode)
    return WellSeparatedIterator(
        testtree, trialtree, trialnode; iswellseparated=f.iswellseparated
    )
end

function WellSeparatedIterator(tree, ::isBlockTree; iswellseparated=iswellseparated)
    return _WellSeparatedIteratorBlockTreeFunctor(iswellseparated)
end

"""
    WellSeparatedIterator(tree, node::Int; iswellseparated=iswellseparated)

Constructs an iterator to identify which translations should occur and which should not.
This determination is based on the concept of well-separated nodes.
Two nodes are considered well-separated if their parents are near each other and
the nodes themselves are far apart.
This assumes that child clusters are completely inside their parent clusters.

# Arguments

  - `tree`: the tree in which the translations occur
  - `node`: the node for which the translations happen
  - `iswellseparated`: a function that returns `true` if two nodes are well-separated and `false` otherwise

# Returns

An iterator that yields the nodes that are well-separated from the specified `node` in the `tree`.
"""
function WellSeparatedIterator(tree, node::Int; iswellseparated=iswellseparated)
    return NodeFilterIterator(tree, node, iswellseparated)
end

"""
    WellSeparatedIterator(testtree, trialtree, trialnode::Int; iswellseparated=iswellseparated)

Constructs an iterator to identify which translations should occur and which should not.
This determination is based on the concept of well-separated nodes.
Two nodes are considered well-separated if their parents are near each other and
the nodes themselves are far apart.
This assumes that child clusters are completely inside their parent clusters.

# Arguments

  - `testtree`: the test tree
  - `trialtree`: the trial tree
  - `trialnode`: the node in the trial tree for which the translations happen
  - `iswellseparated`: a function that returns `true` if two nodes are well-separated and `false` otherwise

# Returns

An iterator that yields the nodes in the `testtree` that are well-separated from the specified `trialnode` in the `trialtree`.
"""
function WellSeparatedIterator(
    testtree, trialtree, trialnode::Int; iswellseparated=iswellseparated
)
    return NodeFilterIterator(testtree, trialtree, trialnode, iswellseparated)
end

function NotWellSeparatedIterator(tree, node::Int; isnotwellseparated=isnotwellseparated)
    return NodeFilterIterator(tree, node, isnotwellseparated)
end

function NotWellSeparatedIterator(
    testtree, trialtree, trialnode::Int; isnotwellseparated=isnotwellseparated
)
    return NodeFilterIterator(testtree, trialtree, trialnode, isnotwellseparated)
end
# Well separated filters ###################################################################

struct _IsWellSeparatedFunctor{IN}
    isnear::IN
end

function (f::_IsWellSeparatedFunctor)(tree)
    return iswellseparated(tree, H2Trees.treetrait(tree); isnear=f.isnear(tree))
end

"""
    iswellseparated
"""
function iswellseparated(; isnear=nothing)
    return _IsWellSeparatedFunctor(isnear)
end

struct _IsWellSeparatedNotBlockTreeFunctor{IN}
    isnear::IN
end

function (f::_IsWellSeparatedNotBlockTreeFunctor)(tree, testnode, trialnode)
    return iswellseparated(
        tree, testnode, trialnode, H2Trees.treetrait(tree); isnear=f.isnear
    )
end

function iswellseparated(tree, ::Any; isnear=isnear)
    return _IsWellSeparatedNotBlockTreeFunctor(isnear)
end

struct _IsWellSeparatedBlockTreeFunctor{IN}
    isnear::IN
end

function (f::_IsWellSeparatedBlockTreeFunctor)(testtree, trialtree, testnode, trialnode)
    return iswellseparated(
        testtree,
        trialtree,
        testnode,
        trialnode,
        treetrait(testtree),
        treetrait(trialtree);
        isnear=f.isnear,
    )
end

function iswellseparated(tree, ::isBlockTree; isnear=isnear)
    return _IsWellSeparatedBlockTreeFunctor(isnear)
end

function iswellseparated(tree, testnode::Int, trialnode::Int)
    return iswellseparated(tree, testnode::Int, trialnode::Int, treetrait(tree))
end

function isnotwellseparated(tree, testnode::Int, trialnode::Int)
    return !iswellseparated(tree, testnode::Int, trialnode::Int)
end

function iswellseparated(testtree, trialtree, testnode::Int, trialnode::Int)
    return iswellseparated(
        testtree,
        trialtree,
        testnode::Int,
        trialnode::Int,
        treetrait(testtree),
        treetrait(trialtree),
    )
end

function isnotwellseparated(testtree, trialtree, testnode::Int, trialnode::Int)
    return !iswellseparated(
        testtree,
        trialtree,
        testnode::Int,
        trialnode::Int,
        treetrait(testtree),
        treetrait(trialtree),
    )
end

function iswellseparated(
    tree, testnode::Int, trialnode::Int, ::AbstractTreeTrait; isnear=isnear
)
    isnear(tree, testnode, trialnode) && return false

    trialparent = H2Trees.parent(tree, trialnode)
    testparent = H2Trees.parent(tree, testnode)

    !isnear(tree, testparent, trialparent) && return false

    return true
end

function iswellseparated(
    testtree,
    trialtree,
    testnode::Int,
    trialnode::Int,
    ::AbstractTreeTrait,
    ::AbstractTreeTrait;
    isnear=isnear,
)
    isnear(testtree, trialtree, testnode, trialnode) && return false

    H2Trees.level(trialtree, trialnode) == H2Trees.minimumlevel(trialtree) && return true
    H2Trees.level(testtree, testnode) == H2Trees.minimumlevel(testtree) && return true

    trialparent = H2Trees.parent(trialtree, trialnode)
    testparent = H2Trees.parent(testtree, testnode)

    !isnear(testtree, trialtree, testparent, trialparent) && return false

    return true
end

# Util functions ###########################################################################

struct _IsTranslatingNodeFunctor{IN}
    TranslatingNodesIterator::IN
end

function (f::_IsTranslatingNodeFunctor)(tree)
    return istranslatingnode(
        tree,
        H2Trees.treetrait(tree);
        TranslatingNodesIterator=f.TranslatingNodesIterator(tree),
    )
end

function istranslatingnode(; TranslatingNodesIterator=nothing)
    return _IsTranslatingNodeFunctor(TranslatingNodesIterator)
end

struct _IsTranslatingNodeNotBlockTreeFunctor{T,IN}
    tree::T
    TranslatingNodesIterator::IN
end

function (f::_IsTranslatingNodeNotBlockTreeFunctor)(node)
    return istranslatingnode(
        f.tree, node; TranslatingNodesIterator=f.TranslatingNodesIterator
    )
end

struct _IsTranslatingNodeBlockTreeFunctor{IN}
    TranslatingNodesIterator::IN
end

function (f::_IsTranslatingNodeBlockTreeFunctor)(testtree, trialtree, trialnode)
    return istranslatingnode(
        testtree, trialtree, trialnode; TranslatingNodesIterator=f.TranslatingNodesIterator
    )
end

function istranslatingnode(tree, ::Any; TranslatingNodesIterator=TranslatingNodesIterator)
    return _IsTranslatingNodeNotBlockTreeFunctor(tree, TranslatingNodesIterator)
end

function istranslatingnode(
    tree, ::isBlockTree; TranslatingNodesIterator=TranslatingNodesIterator
)
    return _IsTranslatingNodeBlockTreeFunctor(TranslatingNodesIterator)
end

function istranslatingnode(
    tree, node::Int; TranslatingNodesIterator=TranslatingNodesIterator
)
    iszero(node) && return false

    for _ in TranslatingNodesIterator(tree, node)
        return true
    end

    return false
end

function istranslatingnode(
    testtree, trialtree, trialnode::Int; TranslatingNodesIterator=TranslatingNodesIterator
)
    iszero(trialnode) && return false

    for _ in TranslatingNodesIterator(testtree, trialtree, trialnode)
        return true
    end

    return false
end

function mintranslationlevel(tree; TranslatingNodesIterator=TranslatingNodesIterator)
    return mintranslationlevel(
        tree, treetrait(tree); TranslatingNodesIterator=TranslatingNodesIterator
    )
end

function mintranslationlevel(
    tree, ::AbstractTreeTrait; TranslatingNodesIterator=TranslatingNodesIterator
)
    for treelevel in levels(tree)
        for node in LevelIterator(tree, treelevel)
            if istranslatingnode(
                tree, node; TranslatingNodesIterator=TranslatingNodesIterator
            )
                return treelevel
            end
        end
    end

    return levels(tree)[end]
end

"""
    TranslatingNodesIterator

    This is a wrapper for the `WellSeparatedIterator`.
"""
TranslatingNodesIterator = WellSeparatedIterator
NotTranslatingNodesIterator = NotWellSeparatedIterator
