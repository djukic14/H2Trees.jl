"""
    BoundingBallTree{N,D,T} <: H2ClusterTree

A cluster tree where nodes are bounded by spheres (balls).

This tree structure uses spherical bounding volumes to organize spatial data hierarchically.
Each node is bounded by a ball with a center and radius.

# Type Parameters

  - `N`: The ambient dimension (coordinate space dimension).
  - `D`: The type of nodes.
  - `T`: The type of the radius.

# Fields

  - `nodes::Vector{Node{D}}`: Vector of nodes comprising the tree.
  - `root::Int`: Index of the root node.
  - `center::SVector{N,T}`: Center of the bounding ball of the tree.
  - `radius::T`: Radius of the bounding ball of the tree.
  - `nodesatlevel::Vector{Vector{Int}}`: Vector of vectors, where each inner vector contains the indices of nodes at a specific level.
"""
struct BoundingBallTree{N,D,T} <: H2ClusterTree
    nodes::Vector{Node{D}}
    root::Int
    center::SVector{N,T}
    radius::T
    nodesatlevel::Vector{Vector{Int}}
end

function BoundingBallTree(
    center, radius; minlevel::Int=1, root::Int=1, balldata=BoundingBallData
)
    rootnode = Node(balldata(Int[], center, radius, minlevel), 0, 0, 0)
    return BoundingBallTree([rootnode], root, center, radius, [Int[]])
end

function Base.eltype(::Union{BoundingBallTree{N,D,T},TwoNTree{N,D,T}}) where {N,D,T}
    return SVector{N,T}
end

H2Trees.treetrait(::Type{BoundingBallTree{N,D,T}}) where {N,D,T} = isBoundingBallTree()
