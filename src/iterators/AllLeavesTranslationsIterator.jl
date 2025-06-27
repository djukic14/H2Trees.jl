struct _IsLeafFunctor{T}
    tree::T
end

function (f::_IsLeafFunctor)(tree, nodea, nodeb)
    return isleaf(tree, nodea)
end

#TODO: add BlockTree case
function AllLeavesTranslationsIterator(tree)
    return TranslatingNodesIterator(; iswellseparated=_IsLeafFunctor)(tree)
end
