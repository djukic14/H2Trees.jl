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
