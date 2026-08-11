# Deterministic exhaustive support enumeration.
# In dimension N the optimal ball has a support of at most N+1 balls.

# Canonicalization.

# Geometric key: center coordinates, then radius.
@inline _key(b::Ball{1,T}) where {T} = (b.center[1], b.radius)
@inline _key(b::Ball{2,T}) where {T} = (b.center[1], b.center[2], b.radius)
@inline _key(b::Ball{3,T}) where {T} = (b.center[1], b.center[2], b.center[3], b.radius)

function _canonicalize(balls::AbstractVector{Ball{N,T}}) where {N,T}
    n = length(balls)
    # `sortperm` is stable, so geometrically identical balls keep original relative order and
    # ties ultimately fall back to original index.
    perm = sortperm(1:n; by=i -> _key(balls[i]))
    sorted = Vector{Ball{N,T}}(undef, n)
    @inbounds for k in 1:n
        sorted[k] = balls[perm[k]]
    end
    return sorted, perm
end

# Candidate acceptance.

# Maximum containment residual, or `nothing` if the candidate misses an input beyond tolerance.
function _global_residual(
    ball::Ball{N,T}, balls::AbstractVector{Ball{N,T}}, tol::Tolerance{T}, mult::T
) where {N,T}
    maxres = typemin(T)
    @inbounds for b in balls
        res = containment_residual(ball, b)
        res > mult * tolerance(tol, _containment_scale(ball, b)) && return nothing
        maxres = max(maxres, res)
    end
    return maxres
end

# Enforce global containment, inflate by any positive residual, and keep the best candidate.
function _consider(
    best::Union{Nothing,SEBBResult{N,T}},
    cand::Union{Nothing,Tuple{SEBBResult{N,T},T}},
    balls::AbstractVector{Ball{N,T}},
    tol::Tolerance{T},
) where {N,T}
    cand === nothing && return best
    result, mult = cand
    maxres = _global_residual(result.ball, balls, tol, mult)
    maxres === nothing && return best

    inflation = max(maxres, zero(T))
    safe_ball = Ball(center(result.ball), radius(result.ball) + inflation)
    safe = _replace_ball(result, safe_ball)
    return _is_better(safe, best, tol) ? safe : best
end

# Deterministic comparison: smallest safe radius, then lexicographic center, then support.
function _is_better(
    cand::SEBBResult{N,T}, best::Union{Nothing,SEBBResult{N,T}}, tol::Tolerance{T}
) where {N,T}
    best === nothing && return true
    rc = radius(cand.ball)
    rb = radius(best.ball)
    tt = tolerance(tol, max(rc, rb))
    rc < rb - tt && return true
    rc > rb + tt && return false

    cc = center(cand.ball)
    cb = center(best.ball)
    @inbounds for k in 1:N
        cc[k] < cb[k] - tt && return true
        cc[k] > cb[k] + tt && return false
    end
    cand.support_size != best.support_size && return cand.support_size < best.support_size
    return cand.support < best.support
end

# Enumeration.

function _enumerate(balls::AbstractVector{Ball{N,T}}, tol::Tolerance{T}) where {N,T}
    n = length(balls)
    best = nothing

    for i in 1:n
        best = _consider(best, _candidate_one(balls, i, tol), balls, tol)
    end

    if n >= 2
        for i in 1:(n - 1), j in (i + 1):n
            best = _consider(best, _candidate_two(balls, i, j, tol), balls, tol)
        end
    end

    if N >= 2 && n >= 3
        for i in 1:(n - 2), j in (i + 1):(n - 1), k in (j + 1):n
            best = _consider(best, _candidate_many(balls, (i, j, k), tol), balls, tol)
        end
    end

    if N >= 3 && n >= 4
        for i in 1:(n - 3), j in (i + 1):(n - 2), k in (j + 1):(n - 1), l in (k + 1):n
            best = _consider(best, _candidate_many(balls, (i, j, k, l), tol), balls, tol)
        end
    end

    return best
end

# Top-level API.

"""
    _smallest_enclosing_ball_with_certificate(balls; atol, rtol) -> SEBBResult

Internal solver returning the enclosing ball together with its optimality certificate. The
support indices in the returned certificate refer to the caller's original input ordering.

Never fails on a valid input: a configuration the exact enumeration cannot resolve falls back
to [`_approximateresult`](@ref), whose certificate is empty. Use [`isexactresult`](@ref) to
tell the two apart.
"""
function _smallest_enclosing_ball_with_certificate(
    balls::AbstractVector{Ball{N,T}};
    atol::Real=_default_atol(T),
    rtol::Real=_default_rtol(T),
) where {N,T}
    isempty(balls) && throw(ArgumentError("at least one ball is required"))
    tol = Tolerance{T}(atol, rtol)

    sorted, perm = _canonicalize(balls)
    best = _enumerate(sorted, tol)
    best === nothing && return _approximateresult(
        balls, "no support set produced a numerically valid tangent ball"
    )
    result = best::SEBBResult{N,T}

    # Remap support indices from canonical order back to the caller's original order.
    orig_support = ntuple(Val(4)) do k
        return k <= result.support_size ? perm[result.support[k]] : 0
    end
    remapped = SEBBResult{N,T}(
        result.ball, orig_support, result.weights, result.support_size
    )

    _validate_result(remapped, balls, tol) ||
        return _approximateresult(balls, "the exact candidate failed its containment check")
    return remapped
end

"""
    _approximateresult(balls, reason) -> SEBBResult

Warn, then return [`_fallbackball`](@ref) wrapped in a certificate-free [`SEBBResult`](@ref).

Failing here instead would abort whatever is being built on top of the solver (a whole tree,
and with it a whole matrix) over a configuration that is merely unresolvable in the caller's
arithmetic. A ball that is guaranteed to enclose but only approximately minimal keeps every
consumer correct, at the price of a marginally larger near field.

The warning uses `maxlog = 1` on purpose: a tree build calls this once per internal node, and a
systematically degenerate input would otherwise emit thousands of identical messages.
"""
function _approximateresult(
    balls::AbstractVector{Ball{N,T}}, reason::AbstractString
) where {N,T}
    @warn "SEBB could not solve this configuration exactly ($reason) and fell back to an approximate enclosing ball. The result is guaranteed to contain every input ball but is not certified minimal. This means the input is numerically degenerate at $T precision; consider a wider float type if the exact ball matters." maxlog =
        1
    z = zero(T)
    return SEBBResult{N,T}(_fallbackball(balls), (0, 0, 0, 0), (z, z, z, z), 0x00)
end

"""
    smallest_enclosing_ball(balls::AbstractVector{Ball{N,T}}; atol, rtol) -> Ball{N,T}

Compute the unique smallest ball enclosing all input balls in ambient dimension 1, 2, or 3.

The algorithm deterministically enumerates support sets of at most `N + 1` balls and solves
the associated tangency equations (Fischer & Gärtner 2004). Negative radii and non-finite
input are rejected at construction time; zero radii are accepted as degenerate point balls.
The returned ball is guaranteed to enclose every input under the arithmetic of the caller
(its radius is inflated by any tiny positive containment residual).

Never raises for a valid nonempty input. A configuration whose tangency algebra cannot be
resolved at `T` precision warns once and returns an approximate (still strictly enclosing)
ball instead; see [`_approximateresult`](@ref).

Tolerances are scale-aware and relative; see [`Tolerance`](@ref).
"""
function smallest_enclosing_ball(
    balls::AbstractVector{Ball{N,T}};
    atol::Real=_default_atol(T),
    rtol::Real=_default_rtol(T),
) where {N,T}
    return _smallest_enclosing_ball_with_certificate(balls; atol=atol, rtol=rtol).ball
end

"""
    smallest_enclosing_ball(centers, radii; atol, rtol) -> Ball

Convenience method taking parallel `centers` and `radii` collections. The ambient dimension
and a common floating-point type are determined once, and a concretely-typed
`Vector{Ball{N,T}}` is built before delegating to the core method.
"""
function smallest_enclosing_ball(centers::AbstractVector, radii::AbstractVector; kwargs...)
    length(centers) == length(radii) ||
        throw(DimensionMismatch("centers and radii must have equal length"))
    isempty(centers) && throw(ArgumentError("at least one ball is required"))

    N = length(first(centers))
    T = float(promote_type(_eltype_of(centers), eltype(radii)))
    balls = Vector{Ball{N,T}}(undef, length(centers))
    @inbounds for i in eachindex(centers)
        balls[i] = Ball(SVector{N,T}(centers[i]), T(radii[i]))
    end
    return smallest_enclosing_ball(balls; kwargs...)
end

# Element type of the coordinate scalars across a centers collection.
_eltype_of(centers) = mapreduce(eltype, promote_type, centers)
