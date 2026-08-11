using H2Trees
using PlotlyJS

const HO = H2Trees.HilbertOrdering

levels = 1:4

# In one dimension the Hilbert curve degenerates to the natural interval order, so the index and
# the coordinate coincide: the traversal is simply left to right at every refinement.
curves = map(levels) do level
    return [HO.hilbertcoordinates(Val(1), i, level)[1] for i in 0:(2 ^ level - 1)]
end

# --- hide-from-docs ---
traces = GenericTrace[]
for (i, level) in enumerate(levels)
    x = curves[i]
    # Normalized to [0, 1] so every level spans the same width and the refinement is comparable
    # row to row.
    push!(
        traces,
        scatter(;
            x=(x .+ 0.5) ./ 2^level,
            y=fill(level, length(x)),
            mode="lines+markers",
            line=attr(; color="rgb(170,170,170)", width=1.5),
            marker=attr(;
                size=max(4, 9 - level), color=0:(length(x) - 1), colorscale="Viridis"
            ),
            name="level $level",
            hovertemplate="index %{marker.color}, cell %{customdata}<extra>level $level</extra>",
            customdata=x,
        ),
    )
end

p = PlotlyJS.plot(
    traces,
    Layout(;
        showlegend=false,
        margin=attr(; l=60, r=20, t=20, b=30),
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis=attr(; title="position along the interval", zeroline=false),
        yaxis=attr(;
            title="",
            tickmode="array",
            tickvals=collect(levels),
            ticktext=["level $level" for level in levels],
            zeroline=false,
        ),
    ),
)
