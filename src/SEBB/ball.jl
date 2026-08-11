# Geometric ball type for the smallest-enclosing-ball-of-balls solver.
#
# This type is intentionally free of any H2Trees tree metadata (no node ids, levels,
# stored values, ...). Keeping it self-contained is what allows the whole `SEBB`
# submodule to be extracted into a standalone package later. Do not couple it to
# `BoundingBallData`.

"""
    Ball{N,T<:AbstractFloat}

An immutable geometric ball in ambient dimension `N` with an `SVector{N,T}` center and a
scalar `radius::T`.

Only dimensions `N ∈ (1, 2, 3)` are supported. Centers must be finite and radii must be
finite and nonnegative. Zero radii are accepted as degenerate point balls.

# Fields

  - `center::SVector{N,T}`
  - `radius::T`
"""
struct Ball{N,T<:AbstractFloat}
    center::SVector{N,T}
    radius::T

    function Ball(center::SVector{N,T}, radius::T) where {N,T<:AbstractFloat}
        N in 1:3 ||
            throw(ArgumentError("SEBB supports ambient dimensions 1, 2, and 3; got $N"))
        all(isfinite, center) ||
            throw(ArgumentError("ball center coordinates must be finite; got $center"))
        isfinite(radius) || throw(ArgumentError("ball radius must be finite; got $radius"))
        radius >= zero(T) ||
            throw(ArgumentError("ball radius must be nonnegative; got $radius"))
        return new{N,T}(center, radius)
    end
end

# Promote a floating center with a possibly non-`T` (but real) radius.
function Ball(center::SVector{N,T}, radius::Real) where {N,T<:AbstractFloat}
    return Ball(center, T(radius))
end

"""
    Ball(center, radius)

Construct a [`Ball`](@ref) from any `AbstractVector` center and real `radius`, promoting the
element type and radius to a common floating-point type.
"""
function Ball(center::AbstractVector, radius::Real)
    N = length(center)
    T = float(promote_type(eltype(center), typeof(radius)))
    return Ball(SVector{N,T}(center), T(radius))
end

# Explicitly-parameterized construction used by the concretely-typed convenience wrapper.
function Ball{N,T}(center, radius) where {N,T<:AbstractFloat}
    return Ball(SVector{N,T}(center), T(radius))
end

"""
    center(ball::Ball)

Return the center of `ball` as an `SVector`.
"""
center(ball::Ball) = ball.center

"""
    radius(ball::Ball)

Return the radius of `ball`.
"""
radius(ball::Ball) = ball.radius

Base.eltype(::Type{Ball{N,T}}) where {N,T} = T
Base.eltype(ball::Ball) = eltype(typeof(ball))

dimension(::Type{Ball{N,T}}) where {N,T} = N
dimension(ball::Ball) = dimension(typeof(ball))

function Base.:(==)(a::Ball, b::Ball)
    return a.center == b.center && a.radius == b.radius
end

function Base.show(io::IO, ball::Ball{N,T}) where {N,T}
    return print(io, "Ball{", N, ",", T, "}(", ball.center, ", ", ball.radius, ")")
end
