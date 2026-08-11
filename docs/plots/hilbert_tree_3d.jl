using H2Trees
using PlotlyJS
using StaticArrays

# A complete 4x4x4 grid, one point per finest cell, so every box on every level is occupied and
# the traversal below is the full Hilbert curve rather than a subsequence of it.
side = 4
points = [
    SVector(i + 0.5, j + 0.5, k + 0.5) for i in 0:(side - 1) for j in 0:(side - 1) for
    k in 0:(side - 1)
]

tree = buildtree(points; builder=TwoNTreeBuilder(; minhalfsize=0.5, minvalues=0))

# `nodesatlevel` is Hilbert-ordered within each level, so node centers trace the curve directly.
levelcurves = map(H2Trees.levels(tree)) do level
    nodes = H2Trees.nodesatlevel(tree, level)
    centers = [H2Trees.center(tree, node) for node in nodes]
    return (getindex.(centers, 1), getindex.(centers, 2), getindex.(centers, 3), nodes)
end

# --- hide-from-docs ---
treelevels = collect(H2Trees.levels(tree))

traces = GenericTrace[]
for (i, level) in enumerate(treelevels)
    x, y, z, ids = levelcurves[i]
    push!(
        traces,
        scatter3d(;
            x=x,
            y=y,
            z=z,
            mode="lines+markers",
            line=attr(; color=0:(length(x) - 1), colorscale="Viridis", width=5),
            marker=attr(; size=4, color=0:(length(x) - 1), colorscale="Viridis"),
            customdata=ids,
            hovertemplate="node %{customdata}<br>(%{x}, %{y}, %{z})<extra></extra>",
            visible=(i == 1),
        ),
    )
end

levellabel(level) = "level $level — $(length(H2Trees.nodesatlevel(tree, level))) nodes"

function titleannotation(level)
    return [
        attr(;
            text=levellabel(level),
            showarrow=false,
            xref="paper",
            yref="paper",
            x=0.5,
            y=1.0,
        ),
    ]
end

steps = [
    attr(;
        label="$level",
        method="update",
        args=[
            # One trace per level here, so the mask is just "this level's trace".
            attr(; visible=[l == level for l in treelevels]),
            attr(; annotations=titleannotation(level)),
        ],
    ) for level in treelevels
]

p = PlotlyJS.plot(
    traces,
    Layout(;
        showlegend=false,
        margin=attr(; l=0, r=0, t=30, b=0),
        # A fixed cubic aspect and a shared camera keep the box in place while stepping levels,
        # so the eye compares paths instead of re-orienting each time.
        scene=attr(;
            xaxis=attr(; visible=false),
            yaxis=attr(; visible=false),
            zaxis=attr(; visible=false),
            aspectmode="cube",
        ),
        annotations=titleannotation(first(treelevels)),
        sliders=[
            attr(;
                active=0,
                currentvalue=attr(; prefix="level: "),
                pad=attr(; t=30),
                steps=steps,
            ),
        ],
    ),
)
