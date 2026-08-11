# Stable real-root solver for the scalar quadratic `A R^2 + B R + C = 0`.
#
# The algebraic roots are *not* a correctness certificate on their own: every returned root
# must still pass the geometric validation in `support.jl`. This solver only guarantees a
# numerically well-behaved set of candidate roots in ascending order.

"""
    QuadraticRoots{T}

Fixed-size result of [`real_roots`](@ref): `count ∈ (0, 1, 2)` valid entries stored in
`roots` in ascending order. Unused entries are zero. A fixed layout avoids a runtime union
over the number of roots.
"""
struct QuadraticRoots{T<:AbstractFloat}
    count::Int
    roots::SVector{2,T}
end

QuadraticRoots{T}() where {T} = QuadraticRoots{T}(0, SVector(zero(T), zero(T)))

# Keep whichever of the two candidate roots is representable. A root can overflow when the
# leading coefficient is a subnormal remnant rather than a real quadratic term; the companion
# root computed through the Citardauq pairing stays accurate in that case.
#
# At most one of them can overflow: their product is `c / a`, and both coefficients are
# normalized to magnitude at most one, so `|r1 * r2| <= 1 / floatmin` (far below what two
# simultaneously-overflowing roots would require). Should that reasoning ever fail, the escaping
# infinity is still caught by `_validate_root`'s own `isfinite` test.
@inline function _finiteroots(r1::T, r2::T) where {T}
    isfinite(r1) && isfinite(r2) && return QuadraticRoots{T}(2, SVector(minmax(r1, r2)...))
    finite = isfinite(r1) ? r1 : r2
    return QuadraticRoots{T}(1, SVector(finite, finite))
end

"""
    real_roots(A, B, C, tol::Tolerance) -> QuadraticRoots

Return the real roots of `A R^2 + B R + C = 0`, robust to degeneracy and cancellation.

The coefficients are normalized by their maximum magnitude. The equation is solved as linear
only when the quadratic coefficient is exactly zero; otherwise the cancellation-resistant
Citardauq pairing is used, which stays accurate for *both* roots however small the quadratic
term becomes, so no threshold has to decide when a quadratic is "really" linear. A slightly
negative discriminant is clamped to zero, coincident roots are deduplicated, and a root that
overflows is dropped. Roots come back in ascending order.
"""
function real_roots(A::T, B::T, C::T, tol::Tolerance{T}) where {T<:AbstractFloat}
    z = zero(T)
    scale = max(abs(A), abs(B), abs(C))
    iszero(scale) && return QuadraticRoots{T}()

    a = A / scale
    b = B / scale
    c = C / scale

    if iszero(a)
        # Genuinely linear: b R + c = 0. The test is exact equality on purpose. `R` is a
        # radius, so `A`, `B` and `C` carry different powers of length; comparing `a` against a
        # relative tolerance would mix units, and did: for a support whose tangent radius is
        # `R` the constant term grows like `R^2` while the quadratic term stays O(1), so a
        # perfectly well-conditioned equation was declared "linear" -- and then rejected, since
        # its linear coefficient is negligible too -- as soon as the geometry grew past roughly
        # `1/sqrt(rtol)` units across. No support of three or more balls then produced a
        # candidate and the whole solve failed.
        iszero(b) && return QuadraticRoots{T}()
        r = -c / b
        return isfinite(r) ? QuadraticRoots{T}(1, SVector(r, r)) : QuadraticRoots{T}()
    end

    disc = b * b - 4 * a * c
    disc < -_unit_tolerance(tol) && return QuadraticRoots{T}()

    s = sqrt(max(disc, z))
    # Cancellation-resistant Citardauq/quadratic pairing.
    q = -(b + copysign(s, b)) / 2
    r1 = q / a
    r2 = iszero(q) ? -b / (2a) : c / q
    roots = _finiteroots(r1, r2)
    roots.count == 2 || return roots

    lo, hi = roots.roots[1], roots.roots[2]
    if abs(hi - lo) <= tolerance(tol, max(abs(lo), abs(hi)))
        mid = (lo + hi) / 2
        return QuadraticRoots{T}(1, SVector(mid, mid))
    end
    return roots
end
