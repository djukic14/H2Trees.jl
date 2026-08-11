using H2Trees
using PlotlyJS

const HO = H2Trees.HilbertOrdering

levels = 1:3

curves = map(levels) do level
    cells = [HO.hilbertcoordinates(Val(3), i, level) for i in 0:(8 ^ level - 1)]
    return (getindex.(cells, 1), getindex.(cells, 2), getindex.(cells, 3))
end

# --- hide-from-docs ---
p = make_subplots(;
    rows=1,
    cols=3,
    specs=[Spec(; kind="scene") Spec(; kind="scene") Spec(; kind="scene")],
    subplot_titles=["level 1 (8 cells)" "level 2 (64 cells)" "level 3 (512 cells)"],
    horizontal_spacing=0.02,
)

for (i, level) in enumerate(levels)
    x, y, z = curves[i]
    add_trace!(
        p,
        scatter3d(;
            x=x,
            y=y,
            z=z,
            mode="lines+markers",
            line=attr(; color=0:(length(x) - 1), colorscale="Viridis", width=4),
            marker=attr(;
                size=max(2, 5 - level), color=0:(length(x) - 1), colorscale="Viridis"
            ),
            hovertemplate="(%{x}, %{y}, %{z})<extra></extra>",
        );
        row=1,
        col=i,
    )
end

# `aspectmode="data"` keeps cells cubic; hiding the axes keeps attention on the path, whose grid
# coordinates are available on hover anyway.
hiddenaxes = attr(;
    xaxis=attr(; visible=false),
    yaxis=attr(; visible=false),
    zaxis=attr(; visible=false),
    aspectmode="data",
)
relayout!(
    p;
    showlegend=false,
    margin=attr(; l=0, r=0, t=30, b=0),
    scene=hiddenaxes,
    scene2=hiddenaxes,
    scene3=hiddenaxes,
)

p
