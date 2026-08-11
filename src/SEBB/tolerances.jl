# Central tolerance policy.
#
# All numerical comparisons in the SEBB core route through this policy. Do NOT scatter raw
# `eps`, `isapprox`, or magic constants elsewhere. Every threshold is derived from a *local
# geometric scale* (radii and relative distances), never from absolute coordinate norms, so
# that the algorithm stays translation invariant.

"""
    Tolerance{T}

A pair of tolerances `(atol, rtol)` used to build scale-aware thresholds through
[`tolerance`](@ref).
"""
struct Tolerance{T<:AbstractFloat}
    atol::T
    rtol::T
end

Tolerance{T}(atol, rtol) where {T} = Tolerance{T}(T(atol), T(rtol))

_default_atol(::Type{T}) where {T<:AbstractFloat} = zero(T)
_default_rtol(::Type{T}) where {T<:AbstractFloat} = sqrt(eps(T))

"""
    tolerance(tol::Tolerance, scale)

Return the scale-aware absolute threshold `atol + rtol * max(scale, 1)`.

`scale` is a *local geometric* quantity (a radius, a relative distance, a normalized
coefficient magnitude), never an absolute coordinate norm. Clamping the scale at one keeps a
sensible floor for tiny/zero-radius configurations.
"""
@inline function tolerance(tol::Tolerance{T}, scale) where {T}
    s = max(T(scale), one(T))
    return tol.atol + tol.rtol * s
end

# Dimensionless threshold for quantities that are naturally O(1) (barycentric weights,
# normalized quadratic coefficients, relative pivots).
@inline _unit_tolerance(tol::Tolerance{T}) where {T} = tolerance(tol, one(T))
