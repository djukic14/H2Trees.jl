"""
    _boundingsphereofpoints(points)

Return the center and radius of the smallest ball that contains all `points`.

Uses Welzl's algorithm with pivoting (`BoundingSphere.WelzlPivot`), which reads the points a few
times whatever their order. `BoundingSphere`'s default, the move-to-front variant, rereads all
preceding points whenever a point lies outside their smallest ball. For points ordered along a
convex curve or surface, such as the vertices of a boundary mesh, almost every point does, and
the default takes O(n²) time.

The pivoting variant treats a point as inside if its squared distance from the center exceeds
the squared radius by at most `eps()`, an absolute tolerance. The points are therefore mapped to
a box of half-size one first, which makes this tolerance relative. Finally, the radius is
widened to the distance of the farthest point, so that no point lies outside the returned ball,
not even by rounding.
"""
function _boundingsphereofpoints(points::AbstractVector{SVector{N,T}}) where {N,T}
    lo = hi = first(points)
    for point in points
        lo = min.(lo, point)
        hi = max.(hi, point)
    end
    offset = (lo + hi) / 2
    halfsize = maximum(hi - lo) / 2
    # Points that all coincide stay at the origin
    scale = iszero(halfsize) ? one(halfsize) : halfsize

    unitpoints = similar(points)
    for i in eachindex(points, unitpoints)
        unitpoints[i] = (points[i] - offset) / scale
    end
    unitcenter, unitradius = BoundingSphere.boundingsphere!(
        unitpoints, BoundingSphere.WelzlPivot()
    )

    center = offset + scale * SVector{N,T}(unitcenter)
    radius = scale * unitradius
    for point in points
        radius = max(radius, norm(point - center))
    end
    return center, radius
end
