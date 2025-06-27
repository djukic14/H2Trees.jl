"""
    BoundingBallTree
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
