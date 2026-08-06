struct _WellSeparatedIteratorBuilder{IN}
    iswellseparated::IN
end

function (builder::_WellSeparatedIteratorBuilder)(tree)
    return WellSeparatedIterator(
        tree, treetrait(tree); iswellseparated=builder.iswellseparated(tree)
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

A callable builder that returns a `WellSeparatedIterator` if provided with a tree.

# Throws

  - `ArgumentError`: if both `isnear` and `iswellseparated` are provided, or if neither is provided.
"""
function WellSeparatedIterator(;
    iswellseparated=nothing, isnear=isnothing(iswellseparated) ? isnear() : nothing
)
    if !((isnear !== nothing) ⊻ (iswellseparated !== nothing))
        throw(ArgumentError("Supply one of (not both) isnear or iswellseparated"))
    end

    filter = if isnothing(iswellseparated)
        H2Trees.iswellseparated(; isnear=isnear)
    else
        iswellseparated
    end

    return _WellSeparatedIteratorBuilder(filter)
end

struct _WellSeparatedSingleTreeBuilder{IN}
    iswellseparated::IN
end

function (builder::_WellSeparatedSingleTreeBuilder)(tree, node)
    return WellSeparatedIterator(tree, node; iswellseparated=builder.iswellseparated)
end

function WellSeparatedIterator(tree, ::AbstractTreeTrait; iswellseparated=iswellseparated)
    return _WellSeparatedSingleTreeBuilder(iswellseparated)
end

struct _WellSeparatedTwoTreeBuilder{IN}
    iswellseparated::IN
end

function (builder::_WellSeparatedTwoTreeBuilder)(iteratedtree, anchortree, anchornode)
    return WellSeparatedIterator(
        iteratedtree, anchortree, anchornode; iswellseparated=builder.iswellseparated
    )
end

function WellSeparatedIterator(tree, ::isBlockTree; iswellseparated=iswellseparated)
    return _WellSeparatedTwoTreeBuilder(iswellseparated)
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

struct _IsWellSeparatedPredicateBuilder{IN}
    isnear::IN
end

function (builder::_IsWellSeparatedPredicateBuilder)(tree)
    return iswellseparated(tree, treetrait(tree); isnear=builder.isnear(tree))
end

"""
    iswellseparated
"""
function iswellseparated(; isnear=nothing)
    return _IsWellSeparatedPredicateBuilder(isnear)
end

struct _IsWellSeparatedSingleTreePredicate{IN}
    isnear::IN
end

function (predicate::_IsWellSeparatedSingleTreePredicate)(tree, testnode, trialnode)
    return iswellseparated(
        tree, testnode, trialnode, treetrait(tree); isnear=predicate.isnear
    )
end

function iswellseparated(tree, ::Any; isnear=isnear)
    return _IsWellSeparatedSingleTreePredicate(isnear)
end

struct _IsWellSeparatedTwoTreePredicate{IN}
    isnear::IN
end

function (predicate::_IsWellSeparatedTwoTreePredicate)(
    iteratedtree, anchortree, iteratednode, anchornode
)
    return iswellseparated(
        iteratedtree,
        anchortree,
        iteratednode,
        anchornode,
        treetrait(iteratedtree),
        treetrait(anchortree);
        isnear=predicate.isnear,
    )
end

function iswellseparated(tree, ::isBlockTree; isnear=isnear)
    return _IsWellSeparatedTwoTreePredicate(isnear)
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

    trialparent = parent(tree, trialnode)
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

    level(trialtree, trialnode) == minimumlevel(trialtree) && return true
    level(testtree, testnode) == minimumlevel(testtree) && return true

    trialparent = parent(trialtree, trialnode)
    testparent = parent(testtree, testnode)

    !isnear(testtree, trialtree, testparent, trialparent) && return false

    return true
end

# Util functions ###########################################################################

struct _IsTranslatingNodePredicateBuilder{IN}
    TranslatingNodesIterator::IN
end

function (builder::_IsTranslatingNodePredicateBuilder)(tree)
    return istranslatingnode(
        tree,
        treetrait(tree);
        TranslatingNodesIterator=builder.TranslatingNodesIterator(tree),
    )
end

function istranslatingnode(; TranslatingNodesIterator=nothing)
    return _IsTranslatingNodePredicateBuilder(TranslatingNodesIterator)
end

struct _IsTranslatingNodeSingleTreePredicate{T,IN}
    tree::T
    TranslatingNodesIterator::IN
end

function (predicate::_IsTranslatingNodeSingleTreePredicate)(node)
    return istranslatingnode(
        predicate.tree, node; TranslatingNodesIterator=predicate.TranslatingNodesIterator
    )
end

struct _IsTranslatingNodeTwoTreePredicate{IN}
    TranslatingNodesIterator::IN
end

function (predicate::_IsTranslatingNodeTwoTreePredicate)(
    iteratedtree, anchortree, anchornode
)
    return istranslatingnode(
        iteratedtree,
        anchortree,
        anchornode;
        TranslatingNodesIterator=predicate.TranslatingNodesIterator,
    )
end

function istranslatingnode(tree, ::Any; TranslatingNodesIterator=TranslatingNodesIterator)
    return _IsTranslatingNodeSingleTreePredicate(tree, TranslatingNodesIterator)
end

function istranslatingnode(
    tree, ::isBlockTree; TranslatingNodesIterator=TranslatingNodesIterator
)
    return _IsTranslatingNodeTwoTreePredicate(TranslatingNodesIterator)
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

"""
    mintranslationlevel(tree; TranslatingNodesIterator=TranslatingNodesIterator)

Return the first tree level that contains at least one translating node.

A node is considered translating if `istranslatingnode(tree, node; TranslatingNodesIterator=...)`
is `true`.

# Returns

The minimum level with a translating node. If no translating node exists,
the last level of `tree` is returned.
"""
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
