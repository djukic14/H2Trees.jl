"""
    QuadPointsTree
"""
struct QuadPointsTree{T} <: H2ClusterTree
    tree::T
end

function uniquepointstreetrait(::QuadPointsTree)
    return NonUniquePoints()
end

@treewrapper QuadPointsTree
