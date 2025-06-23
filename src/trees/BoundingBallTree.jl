struct BoundingBallTree{N,D,T} <: H2ClusterTree
    nodes::Vector{Node{D}}
    root::Int
    center::SVector{N,T}
    radius::T
    nodesatlevel::Vector{Vector{Int}}
end

struct BoundingBallData{N,T}
    values::Vector{Int}
    center::SVector{N,T}
    radius::T
    level::Int
end

function Base.eltype(::Union{BoundingBallTree{N,D,T},TwoNTree{N,D,T}}) where {N,D,T}
    return SVector{N,T}
end

H2Trees.treetrait(::Type{BoundingBallTree{N,D,T}}) where {N,D,T} = isBoundingBallTree()
