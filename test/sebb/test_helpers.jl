# Shared, INDEPENDENT test-side helpers.
#
# Per the plan's Section 19: we must NOT validate SEBB predicates solely with SEBB predicates.
# The residual below is recomputed here from scratch and never calls `SEBB.containment_residual`.

using LinearAlgebra
using StaticArrays

const SB = H2Trees.SEBB

# Independent containment residual: > 0 means `inner` pokes out of `outer`.
function indep_residual(outer, inner)
    return norm(SB.center(inner) .- SB.center(outer)) + SB.radius(inner) - SB.radius(outer)
end

# Assert that `outer` encloses every ball in `balls`, printing a rich message on failure.
function assert_encloses(outer, balls; rtol=nothing, atol=nothing)
    T = typeof(SB.radius(outer))
    rt = rtol === nothing ? sqrt(eps(T)) * 32 : T(rtol)
    at = atol === nothing ? zero(T) : T(atol)
    for inner in balls
        res = indep_residual(outer, inner)
        scale = max(
            SB.radius(outer),
            SB.radius(inner),
            norm(SB.center(inner) .- SB.center(outer)),
            one(T),
        )
        thr = at + rt * scale
        @test res <= thr
        if !(res <= thr)
            @info "containment failure" residual = res scale = scale threshold = thr outer =
                outer inner = inner
        end
    end
end

# Approximate equality of two balls (center + radius) with a scale-aware tolerance.
function balls_approx(a, b; rtol=nothing, atol=nothing)
    T = promote_type(typeof(SB.radius(a)), typeof(SB.radius(b)))
    rt = rtol === nothing ? sqrt(eps(T)) * 64 : T(rtol)
    at = atol === nothing ? zero(T) : T(atol)
    scale = max(SB.radius(a), SB.radius(b), norm(SB.center(a)), norm(SB.center(b)), one(T))
    thr = at + rt * scale
    return abs(SB.radius(a) - SB.radius(b)) <= thr &&
           norm(SB.center(a) .- SB.center(b)) <= thr
end

# Closed-form 1D oracle: union of intervals [c-r, c+r].
function oracle_1d(centers, radii)
    L = minimum(centers[i] - radii[i] for i in eachindex(centers))
    U = maximum(centers[i] + radii[i] for i in eachindex(centers))
    return (L + U) / 2, (U - L) / 2
end

# Construct balls internally tangent to a chosen output ball B(cstar, Rstar): each ball sits at
# cstar - (Rstar - radii[i]) * directions[i], so ‖cstar - center‖ + radius == Rstar exactly (up
# to rounding), independent of anything in the SEBB implementation itself. When the convex hull
# of `directions` contains the origin, Fischer & Gärtner's Lemma 2.2 criterion certifies that
# B(cstar, Rstar) actually *is* the smallest enclosing ball of the constructed balls, giving an
# externally-known expected answer (plan Section 21).
function tangent_balls(cstar, Rstar, directions, radii)
    return [
        SB.Ball(cstar - (Rstar - radii[i]) * directions[i], radii[i]) for
        i in eachindex(directions)
    ]
end

# Regular-tetrahedron unit directions, centered at the origin (their unweighted average is
# exactly zero, so the origin lies in their convex hull). Generic over `T` so the same
# construction can be probed at `Float64` and `BigFloat` precision.
function tetrahedron_dirs(::Type{T}) where {T}
    raw = SVector{3,T}[
        SVector{3,T}(1, 1, 1),
        SVector{3,T}(1, -1, -1),
        SVector{3,T}(-1, 1, -1),
        SVector{3,T}(-1, -1, 1),
    ]
    return [d / norm(d) for d in raw]
end
