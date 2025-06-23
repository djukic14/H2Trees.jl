
"""
    boundingbox(points::Vector{SVector{D, T}})

Returns halfsize and center of bounding box of points. The halfsize is the half of the length of the edge of the bounding box.
"""
function boundingbox(points::AbstractArray{SVector{D,T},1}) where {D,T}
    min_dim = Vector(points[1])
    max_dim = Vector(points[1])

    for i in eachindex(points)
        for j in 1:D
            min_dim[j] = min_dim[j] < points[i][j] ? min_dim[j] : points[i][j]
            max_dim[j] = max_dim[j] > points[i][j] ? max_dim[j] : points[i][j]
        end
    end

    center = MVector{D,T}(undef)

    length_dim = zeros(T, D)
    for j in 1:D
        length_dim[j] = max_dim[j] - min_dim[j]
        center[j] = (max_dim[j] + min_dim[j]) / T(2)
    end

    halflength = maximum(length_dim) / T(2)

    return center, halflength
end
