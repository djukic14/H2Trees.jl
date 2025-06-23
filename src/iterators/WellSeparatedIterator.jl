struct WellSeparatedIteratorFunctor{IN}
    iswellseparated::IN
end

function (f::WellSeparatedIteratorFunctor)(tree)
    return WellSeparatedIterator(
        tree, treetrait(tree); iswellseparated=f.iswellseparated(tree)
    )
end

function WellSeparatedIterator(; isnear=nothing, iswellseparated=nothing)
    if !((isnear !== nothing) ⊻ (iswellseparated !== nothing))
        error("Supply one of (not both) isnear or iswellseparated")
    end

    filter = if isnothing(iswellseparated)
        H2Trees.iswellseparated(; isnear=isnear)
    else
        iswellseparated
    end

    return WellSeparatedIteratorFunctor(filter)
end

struct WellSeparatedIteratorNotBlockTreeFunctor{IN}
    iswellseparated::IN
end

function (f::WellSeparatedIteratorNotBlockTreeFunctor)(tree, node)
    return WellSeparatedIterator(tree, node; iswellseparated=f.iswellseparated)
end

function WellSeparatedIterator(tree, ::Any; iswellseparated=iswellseparated)
    return WellSeparatedIteratorNotBlockTreeFunctor(iswellseparated)
end

struct WellSeparatedIteratorBlockTreeFunctor{IN}
    iswellseparated::IN
end

function (f::WellSeparatedIteratorBlockTreeFunctor)(testtree, trialtree, trialnode)
    return WellSeparatedIterator(
        testtree, trialtree, trialnode; iswellseparated=f.iswellseparated
    )
end

function WellSeparatedIterator(tree, ::isBlockTree; iswellseparated=iswellseparated)
    return WellSeparatedIteratorBlockTreeFunctor(iswellseparated)
end

function WellSeparatedIterator(tree, node::Int; iswellseparated=iswellseparated)
    return NodeFilterIterator(tree, node, iswellseparated)
end

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

struct IsWellSeparatedFunctor{IN}
    isnear::IN
end

function (f::IsWellSeparatedFunctor)(tree)
    return iswellseparated(tree, H2Trees.treetrait(tree); isnear=f.isnear(tree))
end

function iswellseparated(; isnear=nothing)
    return IsWellSeparatedFunctor(isnear)
end

struct IsWellSeparatedNotBlockTreeFunctor{IN}
    isnear::IN
end

function (f::IsWellSeparatedNotBlockTreeFunctor)(tree, testnode, trialnode)
    return iswellseparated(
        tree, testnode, trialnode, H2Trees.treetrait(tree); isnear=f.isnear
    )
end

function iswellseparated(tree, ::Any; isnear=isnear)
    return IsWellSeparatedNotBlockTreeFunctor(isnear)
end

struct IsWellSeparatedBlockTreeFunctor{IN}
    isnear::IN
end

function (f::IsWellSeparatedBlockTreeFunctor)(testtree, trialtree, testnode, trialnode)
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
    return IsWellSeparatedBlockTreeFunctor(isnear)
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

struct IsTranslatingNodeFunctor{IN}
    TranslatingNodesIterator::IN
end

function (f::IsTranslatingNodeFunctor)(tree)
    return istranslatingnode(
        tree,
        H2Trees.treetrait(tree);
        TranslatingNodesIterator=f.TranslatingNodesIterator(tree),
    )
end

function istranslatingnode(; TranslatingNodesIterator=nothing)
    return IsTranslatingNodeFunctor(TranslatingNodesIterator)
end

struct IsTranslatingNodeNotBlockTreeFunctor{T,IN}
    tree::T
    TranslatingNodesIterator::IN
end

function (f::IsTranslatingNodeNotBlockTreeFunctor)(node)
    return istranslatingnode(
        f.tree, node; TranslatingNodesIterator=f.TranslatingNodesIterator
    )
end

struct IsTranslatingNodeBlockTreeFunctor{IN}
    TranslatingNodesIterator::IN
end

function (f::IsTranslatingNodeBlockTreeFunctor)(testtree, trialtree, trialnode)
    return istranslatingnode(
        testtree, trialtree, trialnode; TranslatingNodesIterator=f.TranslatingNodesIterator
    )
end

function istranslatingnode(tree, ::Any; TranslatingNodesIterator=TranslatingNodesIterator)
    return IsTranslatingNodeNotBlockTreeFunctor(tree, TranslatingNodesIterator)
end

function istranslatingnode(
    tree, ::isBlockTree; TranslatingNodesIterator=TranslatingNodesIterator
)
    return IsTranslatingNodeBlockTreeFunctor(TranslatingNodesIterator)
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

TranslatingNodesIterator = WellSeparatedIterator

NotTranslatingNodesIterator = NotWellSeparatedIterator
