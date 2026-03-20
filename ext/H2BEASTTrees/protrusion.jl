"""
    (f::BEASTProtrusionFunctor)(tree, node::Int, value::Int)

Compute the protrusion of basis function `value` relative to `node` in `tree`.

This method uses `H2Trees.center(tree, node)` and `H2Trees.halfsize(tree, node)` and
forwards to the center/halfsize overload.
"""
function (f::BEASTProtrusionFunctor)(tree, node::Int, value::Int)
    return f(H2Trees.center(tree, node), H2Trees.halfsize(tree, node), value)
end

"""
    (f::BEASTProtrusionFunctor)(center::A, halfsize::T, value::Int) where {T,A<:AbstractVector{T}}

Compute the maximum normalized protrusion of basis function `value` with respect to
an axis-aligned box centered at `center` with halfsize `halfsize`.

The protrusion is evaluated over all support vertices and the maximum value is returned.
"""
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
