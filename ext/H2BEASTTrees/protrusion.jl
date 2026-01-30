function (f::BEASTProtrusionFunctor)(tree, node::Int, value::Int)
    return f(H2Trees.center(tree, node), H2Trees.halfsize(tree, node), value)
end

function (f::BEASTProtrusionFunctor)(
    center::A, halfsize::T, value::Int
) where {T,A<:AbstractVector{T}}
    maxprotrusion = zero(T)
    for shape in f.space.fns[value]
        support = BEAST.chart(geometry(f.space), shape.cellid)
        for v in BEAST.vertices(support)
            for i in eachindex(v)
                maxprotrusion = max(
                    maxprotrusion, (abs(v[i] - center[i]) - halfsize) / (2 * halfsize)
                )
            end
        end
    end
    return T(maxprotrusion)
end
