struct QuadPointsTree{T} <: H2ClusterTree
    tree::T
end

function uniquepointstreetrait(tree::QuadPointsTree)
    return NonUniquePoints()
end

@treewrapper QuadPointsTree
