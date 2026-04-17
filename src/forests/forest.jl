"""
        Forest{N,T}

A container for multiple trees.

`Forest` stores a fixed-size collection of trees in an `SVector` and provides
basic indexing and iteration.

# Type Parameters

    - `N`: Number of trees in the forest.
    - `T`: Type of each tree.

# Fields

    - `trees::SVector{N,T}`: Fixed-size collection of trees.
"""
struct Forest{N,T}
    trees::SVector{N,T}
end

function Base.length(f::F) where {F<:Forest}
    return length(f.trees)
end

function Base.getindex(f::F, i) where {F<:Forest}
    return f.trees[i]
end

Base.IteratorSize(::Type{<:Forest}) = Base.HasLength()
Base.eltype(::Type{Forest{N,T}}) where {N,T} = T

function Base.iterate(f::F) where {F<:Forest}
    return iterate(f.trees)
end

function Base.iterate(f::F, state) where {F<:Forest}
    return iterate(f.trees, state)
end
