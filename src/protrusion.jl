function maxprotrusion(tree; computeprotrusion=ComputeProtrusionFunctor)
    T = eltype(eltype(tree))
    # protrusion = Dict(levels(tree) .=> zeros(T, numberoflevels(tree)))
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

# function maxprotrusion(tree; computeprotrusion=ComputeProtrusionFunctor)
#     T = eltype(eltype(tree))
#     protrusion = zeros(T, length(levels(tree)))

#     for level in H2Trees.levels(tree)
#         for node in H2Trees.LevelIterator(tree, level)
#             for val in H2Trees.values(tree, node)
#                 protrusion[H2Trees.leveltolevelid(tree, level)] = max(
#                     protrusion[H2Trees.leveltolevelid(tree, level)],
#                     computeprotrusion(tree, node, val),
#                 )
#             end
#         end
#         protrusion[H2Trees.leveltolevelid(tree, level)] >= 0.5 && @warn(
#             "Max protrusion too large (>= 0.5) at level $level: $(protrusion[H2Trees.leveltolevelid(tree, level)])"
#         )
#     end

#     return protrusion
# end

# relative protrusion normalized to 2*halfsize of the box
struct ComputeProtrusionFunctor end

function (f::ComputeProtrusionFunctor)(tree, leaf::Int, value::Int)
    return f(center(tree, leaf), halfsize(tree, leaf), value)
end

function (f::ComputeProtrusionFunctor)(
    center::A, halfsize::T, value::Int
) where {T,A<:AbstractVector{T}}
    return zero(T)
end

struct MaxProtrusionFunctor{C,M}
    computeprotrusion::C
    maxprotrusion::M
    function MaxProtrusionFunctor(c, m)
        return new{typeof(c),typeof(2 * m)}(c, 2 * m)
    end
end

function MaxProtrusionFunctor(maxprotrusion)
    return MaxProtrusionFunctor(ComputeProtrusionFunctor(), maxprotrusion)
end

function (f::MaxProtrusionFunctor)(center, halfsize, point)
    return f.computeprotrusion(center, halfsize, point) <= f.maxprotrusion
end

# needs H2BEASTTrees extension
struct BEASTProtrusionFunctor{S}
    space::S
end
function MaxBEASTProtrusionFunctor(space, maxprotrusion)
    return MaxProtrusionFunctor(BEASTProtrusionFunctor(space), maxprotrusion)
end
