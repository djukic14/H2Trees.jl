"""
    tracecube(center::SVector{D, T}, halfsize; kwargs...)

Returns a trace, which can be used to plot a bounding box in PlotlyJS. All 'kwargs' are
passed to `PlotlyJS.scatter` or `PlotlyJS.scatter3d`, respectively.

# Arguments:

  - `center`: Center of the bounding box.
  - `halfsize`: Halfsize of the bounding box, which is half of the length of the edge of the bounding box.
  - `kwargs`: keyword arguments passed to `PlotlyJS`
"""
function tracecube(center::SVector{2,T}, halfsize::T; kwargs...) where {T}
    x, y = center
    trace = PlotlyJS.scatter(;
        x=[x - halfsize, x - halfsize, x + halfsize, x + halfsize, x - halfsize],
        y=[y - halfsize, y + halfsize, y + halfsize, y - halfsize, y - halfsize],
        kwargs...,
    )

    return trace
end

function tracecube(center::SVector{3,T}, halfsize::T; kwargs...) where {T}
    x, y, z = center
    trace = PlotlyJS.scatter3d(;
        x=[
            x - halfsize,
            x + halfsize,
            x + halfsize,
            x - halfsize,
            x - halfsize,
            x - halfsize,
            x + halfsize,
            x + halfsize,
            x - halfsize,
            x - halfsize,
            x - halfsize,
            NaN,
            x - halfsize,
            x - halfsize,
            NaN,
            x + halfsize,
            x + halfsize,
            NaN,
            x + halfsize,
            x + halfsize,
        ],
        y=[
            y - halfsize,
            y - halfsize,
            y + halfsize,
            y + halfsize,
            y - halfsize,
            y - halfsize,
            y - halfsize,
            y + halfsize,
            y + halfsize,
            y - halfsize,
            y - halfsize,
            NaN,
            y + halfsize,
            y + halfsize,
            NaN,
            y - halfsize,
            y - halfsize,
            NaN,
            y + halfsize,
            y + halfsize,
        ],
        z=[
            z - halfsize,
            z - halfsize,
            z - halfsize,
            z - halfsize,
            z - halfsize,
            z + halfsize,
            z + halfsize,
            z + halfsize,
            z + halfsize,
            z + halfsize,
            z + halfsize,
            NaN,
            z + halfsize,
            z - halfsize,
            NaN,
            z - halfsize,
            z + halfsize,
            NaN,
            z - halfsize,
            z + halfsize,
        ],
        kwargs...,
    )

    return trace
end

"""
    tracecube(tree::H2ClusterTree, node::Int; kwargs...)

Returns a trace, which can be used to plot a cluster of a `TwoNTree` in `PlotlyJS`. All
'kwargs' are passed to `PlotlyJS.scatter` or `PlotlyJS.scatter3d`, respectively.

# Arguments:

  - `tree::H2ClusterTree`: tree
  - `node::Int`: Cluster to plot
  - `kwargs`: keyword arguments passed to `PlotlyJS`
"""
function tracecube(tree::H2ClusterTree, node::Int; kwargs...)
    return tracecube(center(tree, node), halfsize(tree, node); kwargs...)
end
