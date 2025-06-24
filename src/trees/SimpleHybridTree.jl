# upper tree: [minlevel, hybridlevel]
# lower tree: [hybridlevel+1, maxlevel]
struct SimpleHybridTree{T} <: H2ClusterTree
    tree::T
    hybridlevel::Int
end

function SimpleHybridTree(tree; kwargs...)
    return SimpleHybridTree(tree, treetrait(tree); kwargs...)
end

function SimpleHybridTree(
    tree::T, ::isTwoNTree; hs=H2Trees.halfsizes(tree), hybridhalfsize::H=maximum(hs)
) where {T<:TwoNTree,H<:AbstractFloat}
    hybridhalfsize < minimum(hs) && error(
        "Hybrid halfsize $hybridhalfsize is smaller than minimum halfsize $(minimum(hs))",
    )

    hybridlevel = findfirst(x -> x <= hybridhalfsize, hs)

    isnothing(hybridlevel) ? hybridlevel = H2Trees.root(tree) : hybridlevel -= 1
    hybridlevel = max(hybridlevel, H2Trees.root(tree))

    for leaf in H2Trees.leaves(tree)
        level(tree, leaf) <= hybridlevel &&
            error("Leaf $leaf is below hybrid level $hybridlevel")
    end

    return SimpleHybridTree{T}(tree, hybridlevel)
end

function printtree(io::IO, tree::SimpleHybridTree)
    p = printtree(io, tree.tree)
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
