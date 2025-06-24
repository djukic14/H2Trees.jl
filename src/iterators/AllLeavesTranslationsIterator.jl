struct IsLeafFunctor{T}
    tree::T
end

function (f::IsLeafFunctor)(tree, nodea, nodeb)
    return isleaf(tree, nodea)
end

#TODO: add BlockTree case
function AllLeavesTranslationsIterator(tree)
    return TranslatingNodesIterator(; iswellseparated=IsLeafFunctor)(tree)
end
