using H2Trees
using PlotlyJS

# `HilbertOrdering` is internal, but its index <-> coordinate helpers are exactly what is needed
# to draw the curve itself: walk the level-L index range in order and ask for each cell's
# integer grid coordinates.
const HO = H2Trees.HilbertOrdering

levels = 1:4

curves = map(levels) do level
    cells = [HO.hilbertcoordinates(Val(2), i, level) for i in 0:(4 ^ level - 1)]
    return (first.(cells), last.(cells))
end

# --- hide-from-docs ---
p = make_subplots(;
    rows=2,
    cols=2,
    subplot_titles=[
        "level 1 (4 cells)" "level 2 (16 cells)"
        "level 3 (64 cells)" "level 4 (256 cells)"
    ],
    horizontal_spacing=0.08,
    vertical_spacing=0.12,
)

for (i, level) in enumerate(levels)
    x, y = curves[i]
    row, col = fldmod1(i, 2)
    add_trace!(
        p,
        scatter(;
            x=x,
            y=y,
            mode="lines+markers",
            line=attr(; color="rgb(170,170,170)", width=1.5),
            # Marker colour follows the traversal index, so the direction of travel is readable,
            # not just the shape of the path.
            marker=attr(;
                size=max(3, 7 - level), color=0:(length(x) - 1), colorscale="Viridis"
            ),
            hovertemplate="(%{x}, %{y})<extra></extra>",
        );
        row=row,
        col=col,
    )
end

# Keep cells square; grid coordinates are on hover.
relayout!(
    p;
    showlegend=false,
    margin=attr(; l=20, r=20, t=30, b=20),
    plot_bgcolor="rgba(0,0,0,0)",
    xaxis=attr(; visible=false, scaleanchor="y", scaleratio=1),
    xaxis2=attr(; visible=false, scaleanchor="y2", scaleratio=1),
    xaxis3=attr(; visible=false, scaleanchor="y3", scaleratio=1),
    xaxis4=attr(; visible=false, scaleanchor="y4", scaleratio=1),
    yaxis=attr(; visible=false),
    yaxis2=attr(; visible=false),
    yaxis3=attr(; visible=false),
    yaxis4=attr(; visible=false),
)

p
