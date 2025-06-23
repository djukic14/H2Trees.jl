struct SimpleHybridTree{T} <: H2ClusterTree
    tree::T
    hybridlevel::Int
end

function SimpleHybridTree(tree::T, hybridhalfsize::H) where {T<:TwoNTree,H<:AbstractFloat}
    hybridlevel = findfirst(x -> x < hybridhalfsize, H2Trees.halfsizes(tree))
    isnothing(hybridlevel) ? hybridlevel = H2Trees.root(tree) : hybridlevel -= 1
    hybridlevel = max(hybridlevel, H2Trees.root(tree))

    return SimpleHybridTree{T}(tree, hybridlevel)
end

function printtree(io::IO, tree::SimpleHybridTree)
    p = printtree(io, tree.tree)
    #TODO: color upper and lower levels in other color
    println(io, "upper levels start at: ", tree.hybridlevel)
    return p
end

function hybridlevel(tree::SimpleHybridTree)
    return tree.hybridlevel
end

function uniquepointstreetrait(tree::SimpleHybridTree)
    return uniquepointstreetrait(tree.tree)
end

@treewrapper SimpleHybridTree
