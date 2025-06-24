
"""
    tracecube(center::SVector{D, T}, halfsize; line_color = :black, legendgroup = 0, name = Bounding Box, showlegend = false)

Returns a trace, which can be used to plot a bounding box in PlotlyJS.

# Arguments:

  - `center`: Center of the bounding box.
  - `halfsize`: Halfsize of the bounding box, which is half of the length of the edge of the bounding box.
  - `line_color`: Color of the bounding box in plot.
  - `legendgroup`: Legend group of the trace.
  - `name`: Name of trace. This is used in the legend.
  - `showlegend`: Determines whether or not an item corresponding to this trace is shown in the legend.
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

function tracecube(tree::H2ClusterTree, node::Int; kwargs...)
    return tracecube(center(tree, node), halfsize(tree, node); kwargs...)
end
