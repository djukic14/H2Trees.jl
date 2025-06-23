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

function treetrait(tree::SimpleHybridTree)
    return treetrait(tree.tree)
end

function nodesatlevel(tree::SimpleHybridTree)
    return nodesatlevel(tree.tree)
end

function samelevelnodes(tree::SimpleHybridTree, node::Int)
    return samelevelnodes(tree.tree, node)
end

function root(tree::SimpleHybridTree)
    return root(tree.tree)
end

function (tree::SimpleHybridTree)(node::Int)
    return tree.tree(node)
end

function printtree(io::IO, tree::SimpleHybridTree)
    p = printtree(io, tree.tree)
    #TODO: color upper and lower levels in other color
    println(io, "upper levels start at: ", tree.hybridlevel)
    return p
end

function center(tree::SimpleHybridTree, node::Int)
    return center(tree.tree, node)
end

function halfsize(tree::SimpleHybridTree, node::Int)
    return halfsize(tree.tree, node)
end

function levels(tree::SimpleHybridTree)
    return levels(tree.tree)
end

function leaves(tree::SimpleHybridTree)
    return leaves(tree.tree)
end

function numberoflevels(tree::SimpleHybridTree)
    return numberoflevels(tree.tree)
end

function values(tree::SimpleHybridTree)
    return values(tree.tree)
end

function sector(tree::SimpleHybridTree, node::Int)
    return sector(tree.tree, node)
end

function data(tree::SimpleHybridTree, node::Int)
    return data(tree.tree, node)
end

function parent(tree::SimpleHybridTree, node::Int)
    return parent(tree.tree, node)
end

function nextsibling(tree::SimpleHybridTree, node::Int)
    return nextsibling(tree.tree, node)
end

function firstchild(tree::SimpleHybridTree, node::Int)
    return firstchild(tree.tree, node)
end

function children(tree::SimpleHybridTree, node::Int)
    return children(tree.tree, node)
end

function uniquepointstreetrait(tree::SimpleHybridTree)
    return uniquepointstreetrait(tree.tree)
end

function uniquepointstree(tree::SimpleHybridTree)
    return uniquepointstree(tree.tree)
end

function numberofnodes(tree::SimpleHybridTree)
    return numberofnodes(tree.tree)
end

function Base.eltype(tree::SimpleHybridTree)
    return Base.eltype(tree.tree)
end

function hybridlevel(tree::SimpleHybridTree)
    return tree.hybridlevel
end

function parentcenterminuschildcenter(tree::SimpleHybridTree, node::Int)
    return parentcenterminuschildcenter(tree.tree, node)
end

function oppositesector(tree::SimpleHybridTree, node::Int)
    return oppositesector(tree.tree, node)
end
