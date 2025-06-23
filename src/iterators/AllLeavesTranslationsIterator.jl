struct IsLeafFunctor{T}
    tree::T
end

function (f::IsLeafFunctor)(tree, nodea, nodeb)
    return H2Trees.isleaf(tree, nodea)
end

function AllLeavesTranslationsIterator(tree)
    return H2Trees.TranslatingNodesIterator(; iswellseparated=IsLeafFunctor)(tree)
end
