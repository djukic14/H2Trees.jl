"""
    (f::BEASTProtrusionFunctor)(center::A, halfsize::T, value::Int) where {T,A<:AbstractVector{T}}

Compute the normalized protrusion of basis function `value` from a candidate box.

The box is centered at `center` with scalar `halfsize`. The returned value is the
maximum normalized excess of any support vertex beyond the box boundary,
measured over every coordinate direction. A value below or equal to the
configured `ProtrusionCheck.max` means the basis function fits the candidate
node closely enough for that protrusion policy.
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
