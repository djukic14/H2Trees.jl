"""
    maxprotrusion(tree; computeprotrusion=ComputeProtrusionFunctor)

Compute the maximum protrusion per tree level.

For each leaf, this evaluates `computeprotrusion(tree, leaf, value)` for all values in the
leaf and stores the leaf maximum on its level. The level-wise protrusion is then propagated
upwards so each coarser level is at least half of the next finer level.

A warning is emitted when a level protrusion is greater than or equal to `0.5`.

# Returns

A vector with one protrusion value per level in `levels(tree)`.
"""
function maxprotrusion(tree; computeprotrusion=ComputeProtrusionFunctor)
    T = eltype(eltype(tree))
    protrusion = zeros(T, length(levels(tree)))

    for leaf in H2Trees.leaves(tree)
        leafprotrusion = zero(T)
        for val in H2Trees.values(tree, leaf)
            leafprotrusion = max(leafprotrusion, computeprotrusion(tree, leaf, val))
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

"""
    ComputeProtrusionFunctor

Default protrusion evaluator.

This functor represents relative protrusion normalized by `2 * halfsize` of a box.
The base implementation returns zero and is intended as a lightweight fallback and
extension point.
"""
struct ComputeProtrusionFunctor end

function (f::ComputeProtrusionFunctor)(tree, leaf::Int, value::Int)
    return f(center(tree, leaf), halfsize(tree, leaf), value)
end

function (f::ComputeProtrusionFunctor)(
    center::A, halfsize::T, value::Int
) where {T,A<:AbstractVector{T}}
    return zero(T)
end

# needs H2BEASTTrees extension
"""
    BEASTProtrusionFunctor(space)

Protrusion evaluator backed by a BEAST space.

Concrete call methods for this type are provided by the H2BEASTTrees extension.
"""
struct BEASTProtrusionFunctor{S}
    space::S
end
