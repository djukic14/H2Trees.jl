"""
    maxprotrusion(tree; computeprotrusion=ComputeProtrusionFunctor)

Compute the maximum protrusion per tree level.

For each leaf, this evaluates the protrusion of every stored value relative to
the leaf box and records the maximum on that leaf's level. Specialized
protrusion functors can use the allocation-friendly
`f(center, halfsize, value)` form.

The level-wise values are propagated upward so each coarser level is at least
half of the next finer level. A warning is emitted for values greater than or
equal to `0.5`. The result has one entry per level in `levels(tree)`.
"""
function maxprotrusion(tree; computeprotrusion=ComputeProtrusionFunctor())
    T = eltype(eltype(tree))
    protrusion = zeros(T, length(levels(tree)))

    for leaf in leaves(tree)
        leafprotrusion = zero(T)
        leafcenter = center(tree, leaf)
        leafhalfsize = halfsize(tree, leaf)
        for val in values(data(tree, leaf))
            leafprotrusion = max(
                leafprotrusion, computeprotrusion(leafcenter, leafhalfsize, val)
            )
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
    levelprotrusions(tree; compute=ComputeProtrusionFunctor())

Return the level-wise protrusion vector for `tree`.

This is the public keyword spelling used by tree-building diagnostics; it
delegates to [`maxprotrusion`](@ref).
"""
function levelprotrusions(tree; compute=ComputeProtrusionFunctor())
    return maxprotrusion(tree; computeprotrusion=compute)
end

"""
    protrusionreport(tree; compute=ComputeProtrusionFunctor())

Return the level with the largest protrusion and its value.

The result is a named tuple `(level, value)`.
"""
function protrusionreport(tree; compute=ComputeProtrusionFunctor())
    protrusions = levelprotrusions(tree; compute=compute)
    isempty(protrusions) && return (level=0, value=zero(eltype(eltype(tree))))
    levelid = argmax(protrusions)
    return (level=levels(tree)[levelid], value=protrusions[levelid])
end

"""
    ComputeProtrusionFunctor

Default protrusion evaluator.

The center/halfsize call form represents relative protrusion normalized by
`2 * halfsize` of a box. The base implementation returns zero and is intended
as a lightweight fallback and extension point.
"""
struct ComputeProtrusionFunctor end

"""
    (f::ComputeProtrusionFunctor)(center, halfsize, value)

Return the normalized protrusion of `value` relative to a candidate box.

The default implementation returns zero.
"""
function (f::ComputeProtrusionFunctor)(
    center::A, halfsize::T, value::Int
) where {T,A<:AbstractVector{T}}
    return zero(T)
end

"""
    BEASTProtrusionFunctor(space)

Protrusion evaluator backed by a BEAST space.

Concrete call methods for this type are provided by the H2BEASTTrees extension,
which can evaluate basis-function support against candidate tree boxes.
"""
struct BEASTProtrusionFunctor{S}
    space::S
end
