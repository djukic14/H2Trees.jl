module TestBoundingSphereOfPoints

using H2Trees
using StaticArrays
using LinearAlgebra
using Test

# Points in order along a circle (2D) or a spiral over a sphere (3D) of radius `scale` around
# `offset`, as the vertices of a boundary mesh are ordered. The smallest enclosing ball is the
# circle or sphere itself.
function orderedcircle(n, scale, offset)
    return [offset + scale * SVector(cospi(2i / n), sinpi(2i / n)) for i in 1:n]
end
function orderedsphere(n, scale, offset)
    return map(1:n) do i
        z = 1 - (2i - 1) / n
        ρ = sqrt(1 - z^2)
        φ = i * π * (3 - sqrt(5)) # golden angle
        return offset + scale * SVector(ρ * cos(φ), ρ * sin(φ), z)
    end
end

@testset "Smallest enclosing ball of points" begin
    # Small balls away from the origin: `BoundingSphere`'s algorithms leave points outside some of
    # these balls, by rounding or by the absolute tolerance of the pivoting variant.
    for scale in (1.0, 1e-3, 1e-6), n in (64, 2000)
        for points in (
            orderedcircle(n, scale, SVector(2.0, 3.0)),
            orderedsphere(n, scale, SVector(2.0, 3.0, 4.0)),
        )
            center, radius = H2Trees._boundingsphereofpoints(points)
            @test all(point -> norm(point - center) <= radius, points)
            @test radius ≈ scale rtol = 1e-8
        end
    end

    # Tiny balls: the absolute tolerance of the pivoting variant exceeds their squared radius
    for points in (
        orderedcircle(64, 1e-9, SVector(0.0, 0.0)),
        orderedsphere(64, 1e-9, SVector(0.0, 0.0, 0.0)),
    )
        center, radius = H2Trees._boundingsphereofpoints(points)
        @test all(point -> norm(point - center) <= radius, points)
        @test radius ≈ 1e-9 rtol = 1e-8
    end

    # A single point, and points that all coincide
    for points in ([SVector(2.0, 3.0)], fill(SVector(2.0, 3.0, 4.0), 5))
        @test H2Trees._boundingsphereofpoints(points) == (first(points), 0.0)
    end
end

end # module TestBoundingSphereOfPoints
