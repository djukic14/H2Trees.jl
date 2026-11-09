# Point-read checks: how often construction reads its input points, counted by a point vector
# that counts its element reads. Like allocation counts, read counts measure work without a clock,
# but they do not depend on the Julia version, so these checks are not pinned.

using Test
using H2Trees
using StaticArrays

"""
    CountingPoints(points)

Point vector that counts how often its elements are read. Copies and `similar` vectors
share the count, so the reads of working copies count as well.
"""
struct CountingPoints{P} <: AbstractVector{P}
    points::Vector{P}
    reads::Base.RefValue{Int}
end

CountingPoints(points) = CountingPoints(points, Ref(0))
Base.size(p::CountingPoints) = size(p.points)
function Base.getindex(p::CountingPoints, i::Int)
    p.reads[] += 1
    return p.points[i]
end
Base.setindex!(p::CountingPoints, value, i::Int) = (p.points[i] = value)
Base.copy(p::CountingPoints) = CountingPoints(copy(p.points), p.reads)
function Base.similar(p::CountingPoints, ::Type{T}, dims::Dims{1}) where {T}
    return CountingPoints(similar(p.points, T, dims), p.reads)
end

# Points ordered along a circle (2D) and along a spiral over a sphere (3D), as the vertices of a
# boundary mesh are. Almost every point lies outside the smallest ball of the points before it,
# so the move-to-front variant of Welzl's algorithm (`BoundingSphere`'s default) rereads all
# preceding points for almost every point.
orderedcircle(n) = [SVector(cospi(2i / n), sinpi(2i / n)) for i in 1:n]
function orderedsphere(n)
    return map(1:n) do i
        z = 1 - (2i - 1) / n
        ρ = sqrt(1 - z^2)
        φ = i * π * (3 - sqrt(5)) # golden angle
        return SVector(ρ * cos(φ), ρ * sin(φ), z)
    end
end

const READS_N = 4000
# Measured 16.1 (circle) and 14.6 (sphere) reads per point, and at most 37 for rotated, scaled
# and shifted circles and spheres. The move-to-front variant needs over 1000 here, about n/4.
const READS_PER_POINT_BUDGET = 100

@testset "Smallest enclosing ball of ordered points" begin
    for points in (orderedcircle(READS_N), orderedsphere(READS_N))
        counted = CountingPoints(points)
        H2Trees._boundingsphereofpoints(counted)
        @test counted.reads[] <= READS_PER_POINT_BUDGET * length(points)
    end
end

@testset "Root ball of a BoundingBallTree of ordered points" begin
    # With `minvalues` at the number of points nothing is split: the root ball is the only part
    # of the construction that reads the points.
    for points in (orderedcircle(READS_N), orderedsphere(READS_N))
        counted = CountingPoints(points)
        perf_buildboundingballtree(counted; minvalues=length(points))
        @test counted.reads[] <= READS_PER_POINT_BUDGET * length(points)
    end
end
