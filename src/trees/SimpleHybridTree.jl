# upper tree: [minlevel, hybridlevel]
# lower tree: [hybridlevel+1, maxlevel]

"""
    SimpleHybridTree{T} <: H2ClusterTree

A cluster tree that divides into upper and lower tree regions at a hybrid level.

The upper tree spans from the minimum level to `hybridlevel`, while the lower tree spans from
`hybridlevel+1` to the maximum level. This structure enables efficient hybrid algorithms that
treat different tree regions differently.

# Fields

  - `tree::T`: The underlying cluster tree.
  - `hybridlevel::Int`: The level that separates the upper and lower tree regions.
"""
struct SimpleHybridTree{T} <: H2ClusterTree
    tree::T
    hybridlevel::Int
end

function SimpleHybridTree(tree; kwargs...)
    return SimpleHybridTree(tree, treetrait(tree); kwargs...)
end

"""
    SimpleHybridTree(tree::TwoNTree; hs=H2Trees.halfsizes(tree), hybridhalfsize::H=maximum(hs))

Construct a `SimpleHybridTree` for a `TwoNTree` by specifying the hybrid halfsize.

The hybrid level is determined by finding the first level whose halfsize is less than or equal
to the specified `hybridhalfsize`. An error is raised if the hybrid halfsize is smaller than
the minimum halfsize in the tree or if any leaf is below the computed hybrid level.

# Arguments

  - `tree::T`: The underlying `TwoNTree` cluster tree.
  - `hs`: The halfsizes of tree levels (default: `H2Trees.halfsizes(tree)`).
  - `hybridhalfsize::H`: The halfsize threshold for determining the hybrid level (default: maximum halfsize).

# Returns

A new `SimpleHybridTree` instance with the computed hybrid level.

# Throws

Error if `hybridhalfsize` is smaller than the minimum halfsize or if any leaf is below the hybrid level.
"""
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

@treewrapper SimpleHybridTree
