"""
    QuadPointsTree
"""
struct QuadPointsTree{T} <: H2ClusterTree
    tree::T
end

function uniquepointstreetrait(::QuadPointsTree)
    return NonUniquePoints()
end

function numberofvalues(tree::QuadPointsTree)
    maxvalue = 0
    for leaf in H2Trees.leaves(tree)
        maxvalue = max(maxvalue, maximum(H2Trees.values(tree, leaf)))
    end
    return maxvalue
end
@treewrapper QuadPointsTree
