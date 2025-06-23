function H2Trees.traceball(center, radius; n=30, kwargs...)
    return traceball(center, radius; n=n, kwargs...)
end

function traceball(center::SVector{3,T}, radius; n=30, kwargs...) where {T}
    x, y, z = sphere(center, radius; n=n)
    return PlotlyJS.surface(; x=x, y=y, z=z, kwargs...)
end

function sphere(center::SVector{3,T}, radius; n=30) where {T}   # r: radius; C: center [cx,cy,cz]
    u = Base.range(-π, π; length=n)
    v = Base.range(0, π; length=n)
    x = center[1] .+ radius * cos.(u) * sin.(v)'
    y = center[2] .+ radius * sin.(u) * sin.(v)'
    z = center[3] .+ radius * ones(n) * cos.(v)'
    return x, y, z
end

function traceball(center::SVector{2,T}, radius; n=30, kwargs...) where {T}
    x, y = sphere(center, radius; n=n)
    return scatter(; x=x, y=y, kwargs...)
end

function sphere(center::SVector{2,T}, radius; n=30) where {T}   # r: radius; C: center [cx,cy]
    u = Base.range(-π, π; length=n)
    x = center[1] .+ radius * cos.(u)'
    y = center[2] .+ radius * sin.(u)'
    return x, y
end
