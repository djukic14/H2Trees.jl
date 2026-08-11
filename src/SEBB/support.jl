# Candidate balls from support sets of size 1..4.
#
# For a support (a subset of the input balls) we compute the ball internally tangent to all
# support balls whose center lies in the affine hull of their centers, together with the
# barycentric weights of the center. Fischer & Gärtner (2004):
#   * Lemma 2.2: a ball internally tangent to a set V is MB(V) iff its center lies in the
#     convex hull of the centers of V  ->  we require all weights >= 0.
#   * Lemma 2.5: the centers of a basis are affinely independent  ->  we reject
#     rank-deficient supports rather than regularizing them.
#   * Lemma 3.1: the tangent ball of a basis is obtained by a small linear system followed by
#     a quadratic equation  ->  `_candidate_many`.

"""
    SEBBResult{N,T}

Internal result carrying the enclosing `ball` together with an optimality certificate: the
`support` indices (into the caller's original ordering, unused entries zero), the barycentric
`weights` of the center over the support centers (unused entries zero), and `support_size`.

The public API returns only `.ball`; this richer form is used by tests and diagnostics.
"""
struct SEBBResult{N,T<:AbstractFloat}
    ball::Ball{N,T}
    support::NTuple{4,Int}
    weights::NTuple{4,T}
    support_size::UInt8
end

"""
    active_support(result::SEBBResult) -> Vector{Int}

Return the active (nonzero-slot) support indices of `result`.
"""
function active_support(result::SEBBResult)
    return [result.support[k] for k in 1:(result.support_size)]
end

"""
    isexactresult(result::SEBBResult) -> Bool

Whether `result` carries an optimality certificate, i.e. whether it came from the exact
enumeration rather than from the approximate [`_fallbackball`](@ref).

A fallback result has an empty support, which is the only way a caller can tell the two apart:
both are guaranteed to enclose every input, but only an exact result is guaranteed minimal.
"""
isexactresult(result::SEBBResult) = !iszero(result.support_size)

"""
    active_weights(result::SEBBResult) -> Vector

Return the active barycentric weights of `result`.
"""
function active_weights(result::SEBBResult{N,T}) where {N,T}
    return T[result.weights[k] for k in 1:(result.support_size)]
end

# Replace the ball of a result, keeping the certificate (used for safe-radius inflation).
@inline function _replace_ball(r::SEBBResult{N,T}, ball::Ball{N,T}) where {N,T}
    return SEBBResult{N,T}(ball, r.support, r.weights, r.support_size)
end

_pad4(t::NTuple{1,Int}) = (t[1], 0, 0, 0)
_pad4(t::NTuple{2,Int}) = (t[1], t[2], 0, 0)
_pad4(t::NTuple{3,Int}) = (t[1], t[2], t[3], 0)
_pad4(t::NTuple{4,Int}) = t

# Requires a nonempty tuple so `T` is always bound (avoids an unbound type parameter).
@inline function _pad4w(w::Tuple{T,Vararg{T}}) where {T}
    z = zero(T)
    m = length(w)
    return ntuple(k -> k <= m ? w[k] : z, Val(4))
end

# One-ball support.

"""
    _candidate_one(balls, i, tol) -> SEBBResult

The one-ball candidate is the input ball itself (weight one). Global feasibility is decided
later by the enumeration. Handles single balls, nested balls, coincident centers, and
duplicate largest balls.
"""
function _candidate_one(
    balls::AbstractVector{Ball{N,T}}, i::Int, ::Tolerance{T}
) where {N,T}
    result = SEBBResult{N,T}(
        balls[i], (i, 0, 0, 0), (one(T), zero(T), zero(T), zero(T)), 0x01
    )
    # A one-ball candidate is exact; its containment check needs no conditioning slack.
    return (result, one(T))
end

# Two-ball support.

"""
    _candidate_two(balls, i, j, tol) -> Union{Nothing,SEBBResult}

Exact two-ball tangent ball. Rejects coincident centers (a one-ball support handles those)
and nested pairs (a barycentric weight leaves `[0, 1]`).
"""
function _candidate_two(
    balls::AbstractVector{Ball{N,T}}, i::Int, j::Int, tol::Tolerance{T}
) where {N,T}
    c1 = center(balls[i])
    r1 = radius(balls[i])
    c2 = center(balls[j])
    r2 = radius(balls[j])

    diff = c2 - c1
    d = norm(diff)
    # Coincident centers: reject, the one-ball support covers this.
    d <= tolerance(tol, max(r1, r2)) && return nothing

    R = (d + r1 + r2) / 2
    t = (R - r1) / d
    weight_tol = _unit_tolerance(tol)
    λ1 = one(T) - t
    λ2 = t
    (λ1 >= -weight_tol && λ2 >= -weight_tol) || return nothing
    R >= max(r1, r2) - tolerance(tol, max(r1, r2)) || return nothing

    c = c1 + t * diff
    ball = Ball(c, R)

    # Verify tangency to both support balls.
    _tangent_ok(ball, balls[i], tol) || return nothing
    _tangent_ok(ball, balls[j], tol) || return nothing

    result = SEBBResult{N,T}(ball, (i, j, 0, 0), (λ1, λ2, zero(T), zero(T)), 0x02)
    # The two-ball tangent ball is exact (closed form); no conditioning slack needed.
    return (result, one(T))
end

@inline function _tangent_ok(
    ball::Ball{N,T}, support::Ball{N,T}, tol::Tolerance{T}, mult::T=one(T)
) where {N,T}
    residual = containment_residual(ball, support)
    return abs(residual) <= mult * tolerance(tol, _containment_scale(ball, support))
end

# Three/four-ball support (Lemma 3.1).

"""
    _candidate_many(balls, idx::NTuple{M,Int}, tol) -> Union{Nothing,SEBBResult}

Tangent ball of a size-`M` support (`M ∈ (3, 4)`) via the Fischer & Gärtner (2004) Lemma 3.1
construction: a Gram linear system parameterized by the radius, solved with a Cholesky
factorization (never an explicit inverse), followed by a quadratic in the radius.

Returns `nothing` for affinely dependent (rank-deficient) supports and for roots that fail
the radius / convex-weight / tangency checks.
"""
function _candidate_many(
    balls::AbstractVector{Ball{N,T}}, idx::NTuple{M,Int}, tol::Tolerance{T}
) where {N,T,M}
    c1 = center(balls[idx[1]])
    r1 = radius(balls[idx[1]])

    # Direction matrix V (N x (M-1)) and Gram matrix G ((M-1) x (M-1)).
    V = _direction_matrix(balls, idx, c1, Val(M))
    G = transpose(V) * V

    # Affine-independence / conditioning guard on the Gram matrix.
    gscale = _maxdiag(G)
    iszero(gscale) && return nothing
    F = cholesky(Symmetric(G); check=false)
    issuccess(F) || return nothing
    minpiv2 = _minpivotsq(F)
    # Reject only genuinely rank-deficient (affinely dependent) supports; a smaller subset
    # then represents the same basis if one exists. Pivots and `gscale` are squared-length
    # quantities, so the comparison is done in linear-length units to route through the same
    # centralized `tolerance(tol, scale)` policy (atol included) used everywhere else.
    sqrt(minpiv2) > tolerance(tol, sqrt(gscale)) || return nothing

    e = _rhs_e(balls, idx, V, r1, Val(M))
    f = _rhs_f(balls, idx, r1, Val(M))

    a = F \ e
    b = F \ f

    # Quadratic A R^2 + B R + C = 0 from ‖V x‖^2 = (R - r1)^2, x = a + R b.
    Ga = G * a
    Gb = G * b
    A = dot(b, Gb) - one(T)
    B = 2 * (dot(a, Gb) + r1)
    C = dot(a, Ga) - r1 * r1

    roots = real_roots(A, B, C, tol)
    roots.count == 0 && return nothing

    # Conditioning-aware relaxation factor for the internal tangency check. A near-affinely-
    # dependent support (near-collinear/near-coplanar centers) produces a near-double quadratic
    # root, so the tangent radius R can only be resolved to ~sqrt(cond)*eps*scale. The support
    # balls are then *contained* (residual <= 0) but not tangent to machine precision. Global
    # containment plus safe-radius inflation remain the hard correctness guarantee (they are
    # enforced in `_consider`), so here we only need to avoid rejecting a numerically valid
    # candidate whose tangency residual reflects the support's conditioning, not a bad root.
    condmult = sqrt(gscale / minpiv2)

    rmax = _maxradius(balls, idx, Val(M))
    for ri in 1:(roots.count)
        R = roots.roots[ri]
        result = _validate_root(balls, idx, V, c1, r1, a, b, R, rmax, tol, condmult)
        result === nothing && continue
        # Ascending roots: first accepted root has the smallest radius for this support. The
        # conditioning multiplier travels with the candidate so the global containment check
        # (and its safe-radius inflation) applies the same conditioning-aware slack.
        return (result, condmult)
    end
    return nothing
end

# Build V column by column. Unrolled per M for type stability of the static matrix.
@inline function _direction_matrix(
    balls, idx::NTuple{3,Int}, c1::SVector{N,T}, ::Val{3}
) where {N,T}
    v1 = center(balls[idx[2]]) - c1
    v2 = center(balls[idx[3]]) - c1
    return hcat(v1, v2)
end
@inline function _direction_matrix(
    balls, idx::NTuple{4,Int}, c1::SVector{N,T}, ::Val{4}
) where {N,T}
    v1 = center(balls[idx[2]]) - c1
    v2 = center(balls[idx[3]]) - c1
    v3 = center(balls[idx[4]]) - c1
    return hcat(v1, v2, v3)
end

@inline function _rhs_e(balls, idx, V, r1::T, ::Val{M}) where {T,M}
    return SVector{M - 1,T}(
        ntuple(Val(M - 1)) do j
            vj = V[:, j]
            rj = radius(balls[idx[j + 1]])
            return (dot(vj, vj) - rj * rj + r1 * r1) / 2
        end,
    )
end

@inline function _rhs_f(balls, idx, r1::T, ::Val{M}) where {T,M}
    return SVector{M - 1,T}(
        ntuple(Val(M - 1)) do j
            return radius(balls[idx[j + 1]]) - r1
        end,
    )
end

@inline function _maxradius(balls, idx::NTuple{M,Int}, ::Val{M}) where {M}
    r = radius(balls[idx[1]])
    for k in 2:M
        r = max(r, radius(balls[idx[k]]))
    end
    return r
end

@inline function _maxdiag(G)
    m = zero(eltype(G))
    for k in 1:size(G, 1)
        m = max(m, abs(G[k, k]))
    end
    return m
end

# Smallest squared pivot of the Cholesky factor. A tiny value (relative to the Gram scale)
# signals near affine dependence. No `det`, no explicit inverse.
@inline function _minpivotsq(F)
    U = F.U
    T = eltype(U)
    minpiv2 = typemax(T)
    for k in 1:size(U, 1)
        minpiv2 = min(minpiv2, U[k, k]^2)
    end
    return minpiv2
end

function _validate_root(
    balls::AbstractVector{Ball{N,T}},
    idx::NTuple{M,Int},
    V,
    c1::SVector{N,T},
    r1::T,
    a,
    b,
    R::T,
    rmax::T,
    tol::Tolerance{T},
    condmult::T,
) where {N,T,M}
    # A radius is nonnegative by definition and `Ball` enforces that, so the `rmax` test alone
    # is not enough: for a geometry much smaller than one, `tolerance`'s `max(scale, 1)` floor
    # makes `rmax - tolerance(tol, rmax)` negative, and a negative root would then reach the
    # `Ball` constructor below and throw instead of simply being rejected here.
    (isfinite(R) && R >= max(rmax - tolerance(tol, rmax), zero(T))) || return nothing

    x = a + R * b
    c = c1 + V * x
    all(isfinite, c) || return nothing
    ball = Ball(c, R)

    weight_tol = _unit_tolerance(tol)
    sumx = sum(x)
    λ1 = one(T) - sumx
    λ1 >= -weight_tol || return nothing
    for j in 1:(M - 1)
        x[j] >= -weight_tol || return nothing
    end

    # Support tangency (tolerance widened by the support's conditioning; see `_candidate_many`).
    for k in 1:M
        _tangent_ok(ball, balls[idx[k]], tol, condmult) || return nothing
    end

    weights = _pad4w((λ1, ntuple(j -> x[j], Val(M - 1))...))
    return SEBBResult{N,T}(ball, _pad4(idx), weights, UInt8(M))
end
