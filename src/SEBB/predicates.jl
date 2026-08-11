# Containment / tangency predicates.
#
# `containment_residual` is the single source of truth for "how far is `inner` from being
# enclosed by `outer`". A positive residual means `inner` pokes out of `outer`.

"""
    containment_residual(outer::Ball, inner::Ball)

Return `‖c_outer - c_inner‖ + r_inner - r_outer`.

This is nonpositive exactly when `inner` is contained in `outer` (up to rounding). The same
expression is the internal-tangency residual: it is zero when `inner` is internally tangent
to `outer`.
"""
function containment_residual(outer::Ball{N,T}, inner::Ball{N,T}) where {N,T}
    return norm(center(outer) - center(inner)) + radius(inner) - radius(outer)
end

# Scale used for containment / tangency comparisons of a specific (outer, inner) pair.
@inline function _containment_scale(outer::Ball{N,T}, inner::Ball{N,T}) where {N,T}
    return max(radius(outer), radius(inner), norm(center(outer) - center(inner)))
end

"""
    encloses(outer::Ball, inner::Ball; atol, rtol) -> Bool

Return `true` when `inner` is contained in `outer` within the scale-aware tolerance.

The containment residual `‖c_outer - c_inner‖ + r_inner - r_outer` must not exceed
`atol + rtol * max(scale, 1)`, where `scale = max(r_outer, r_inner, ‖c_outer - c_inner‖)`.
The scale is purely relative, preserving translation invariance.
"""
function encloses(
    outer::Ball{N,T},
    inner::Ball{N,T};
    atol::Real=_default_atol(T),
    rtol::Real=_default_rtol(T),
) where {N,T}
    tol = Tolerance{T}(atol, rtol)
    residual = containment_residual(outer, inner)
    return residual <= tolerance(tol, _containment_scale(outer, inner))
end
