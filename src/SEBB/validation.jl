# Independent final validation of a result.
#
# This recomputes containment directly against the *original* input balls (not the canonical
# copy) as a last safety net before returning. It guards against remapping mistakes and any
# candidate that slipped through with a stale ball.

"""
    _validate_result(result::SEBBResult, balls, tol::Tolerance) -> Bool

Return `true` when `result.ball` encloses every ball in `balls` within tolerance.
"""
function _validate_result(
    result::SEBBResult{N,T}, balls::AbstractVector{Ball{N,T}}, tol::Tolerance{T}
) where {N,T}
    ball = result.ball
    isfinite(radius(ball)) || return false
    all(isfinite, center(ball)) || return false
    @inbounds for b in balls
        res = containment_residual(ball, b)
        res > tolerance(tol, _containment_scale(ball, b)) && return false
    end
    return true
end
