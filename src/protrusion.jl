function maxprotrusion(tree; computeprotrusion=H2Trees.computeprotrusion)
    T = eltype(eltype(tree))
    # protrusion = Dict(levels(tree) .=> zeros(T, numberoflevels(tree)))
    protrusion = zeros(T, length(levels(tree)))

    for leaf in H2Trees.leaves(tree)
        leafprotrusion = zero(T)
        for val in H2Trees.values(tree, leaf)
            leafprotrusion = maximum(computeprotrusion(tree, leaf, val))
        end

        levelid = leveltolevelid(tree, level(tree, leaf))
        protrusion[levelid] = max(protrusion[levelid], leafprotrusion)
    end

    for level in reverse(levels(tree))
        level == maximum(levels(tree)) && continue
        levelid = leveltolevelid(tree, level)
        protrusion[levelid] = max(protrusion[levelid], protrusion[levelid + 1] / 2)
    end

    for level in levels(tree)
        levelid = leveltolevelid(tree, level)
        protrusion[levelid] >= 0.5 && @warn(
            "Max protrusion too large (>= 0.5) at level $level: $(protrusion[levelid])"
        )
    end

    return protrusion
end

# relative protrusion normalized to 2*halfsize of the box
function computeprotrusion(tree, leaf::Int, value::Int)
    return computeprotrusion(center(tree, leaf), halfsize(tree, leaf), value)
end

function computeprotrusion(
    center::A, halfsize::T, value::Int
) where {T,A<:AbstractVector{T}}
    return zero(T)
end

# needs H2BEASTTrees extension
struct BEASTProtrusionFunctor{S}
    space::S
end
