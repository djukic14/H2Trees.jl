# Approximate fallback for numerically degenerate inputs.
#
# `_enumerate` is exact by construction, so a configuration it rejects is one whose tangency
# algebra cannot be resolved in the caller's arithmetic at all. Failing there would turn a
# rounding accident into a hard failure of whatever is being built on top: a whole tree, and
# with it a whole matrix. A slightly-too-large ball keeps every consumer correct instead: a
# bounding ball is only ever required to CONTAIN its children, never to be minimal, and the
# only cost of a larger one is a slightly larger near field.
#
# The construction below is therefore chosen for its guarantee rather than for optimality. The
# radius returned is the measured enclosing radius at the returned center, so containment holds
# by definition rather than by tolerance; only minimality is approximate.

"""
    _enclosingradius(center, balls)

The exact radius a ball at `center` needs in order to contain every ball in `balls`, namely
`max(‖center - cᵢ‖ + rᵢ)`.
"""
function _enclosingradius(c::SVector{N,T}, balls::AbstractVector{Ball{N,T}}) where {N,T}
    R = zero(T)
    @inbounds for b in balls
        R = max(R, norm(c - center(b)) + radius(b))
    end
    return R
end

# Number of Bădoiu-Clarkson refinement steps. The scheme's error decays like `1/steps`
# It runs only on the rare degenerate input, and each step is `O(length(balls))` with no
# allocation, so a generous count is cheaper than a second failure mode.
const _FALLBACKSTEPS = 512

# Replace the incumbent when centering at `c` needs a strictly smaller radius.
@inline function _considercenter(
    best::Tuple{SVector{N,T},T}, c::SVector{N,T}, balls::AbstractVector{Ball{N,T}}
) where {N,T}
    R = _enclosingradius(c, balls)
    return R < best[2] ? (c, R) : best
end

@inline function _firstaxis(::Type{SVector{N,T}}) where {N,T}
    return SVector{N,T}(ntuple(d -> d == 1 ? one(T) : zero(T), Val(N)))
end

# Exact smallest ball of the two balls that are farthest apart. For two inputs this is already
# the answer, and for more it is a much better starting point than any single input center.
function _diametralcenter(balls::AbstractVector{Ball{N,T}}) where {N,T}
    n = length(balls)
    besti, bestj, bestspan = 1, 1, typemin(T)
    @inbounds for i in 1:(n - 1), j in (i + 1):n
        span =
            norm(center(balls[i]) - center(balls[j])) + radius(balls[i]) + radius(balls[j])
        if span > bestspan
            besti, bestj, bestspan = i, j, span
        end
    end
    besti == bestj && return center(balls[besti])

    c1, r1 = center(balls[besti]), radius(balls[besti])
    diff = center(balls[bestj]) - c1
    d = norm(diff)
    iszero(d) && return c1
    R = (d + r1 + radius(balls[bestj])) / 2
    # Clamped because a nested pair puts the tangent point outside the segment, in which case
    # the containing ball's own center is the right answer.
    return c1 + clamp((R - r1) / d, zero(T), one(T)) * diff
end

# One Bădoiu-Clarkson step: move a `1/(step + 1)` fraction of the way toward the point of the
# currently farthest ball.
function _towardfarthest(
    c::SVector{N,T}, balls::AbstractVector{Ball{N,T}}, step::Int
) where {N,T}
    far = firstindex(balls)
    fardistance = typemin(T)
    @inbounds for k in eachindex(balls)
        d = norm(c - center(balls[k])) + radius(balls[k])
        if d > fardistance
            fardistance, far = d, k
        end
    end
    # The point of `balls[far]` farthest from `c` lies on the far side of its center, so the
    # direction runs from `c` outwards through that center -- not the other way round.
    diff = center(balls[far]) - c
    d = norm(diff)
    direction = iszero(d) ? _firstaxis(SVector{N,T}) : diff / d
    surface = center(balls[far]) + radius(balls[far]) * direction
    return c + (surface - c) / (step + 1)
end

"""
    _fallbackball(balls) -> Ball

A ball that provably encloses every ball in `balls` and is approximately the smallest such.

Used by [`smallest_enclosing_ball`](@ref) when the exact enumeration cannot resolve the
configuration. Deterministic, allocation free, and independent of the input order beyond
tie-breaking: it seeds from the best input center and the exact two-ball solution of the
farthest-apart pair, then runs `_FALLBACKSTEPS` Bădoiu-Clarkson steps, keeping the best center
any of them produced.

The returned radius is `_enclosingradius` evaluated at the returned center, i.e. exactly the
quantity a containment check recomputes, so containment cannot fail by a rounding margin.
"""
function _fallbackball(balls::AbstractVector{Ball{N,T}}) where {N,T}
    first_center = center(first(balls))
    best = (first_center, _enclosingradius(first_center, balls))
    @inbounds for k in 2:length(balls)
        best = _considercenter(best, center(balls[k]), balls)
    end
    best = _considercenter(best, _diametralcenter(balls), balls)

    c = best[1]
    for step in 1:_FALLBACKSTEPS
        c = _towardfarthest(c, balls, step)
        best = _considercenter(best, c, balls)
    end
    return Ball(best[1], best[2])
end
