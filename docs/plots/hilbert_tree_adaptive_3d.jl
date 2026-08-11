using H2Trees
using PlotlyJS
using StaticArrays
using Random

# Dense core plus sparse halo, as in the 2D adaptive figure: `minvalues` stops subdivision early
# in the halo and late in the core, so leaves end up spread over many levels.
Random.seed!(23)
points = vcat(
    [SVector((0.05 .* randn(3))...) for _ in 1:800], [SVector(randn(3)...) for _ in 1:250]
)

tree = buildtree(points; builder=TwoNTreeBuilder(; minvalues=20))

levelcurves = map(H2Trees.levels(tree)) do level
    nodes = H2Trees.nodesatlevel(tree, level)
    centers = [H2Trees.center(tree, node) for node in nodes]
    return (getindex.(centers, 1), getindex.(centers, 2), getindex.(centers, 3), nodes)
end

# --- hide-from-docs ---
treelevels = collect(H2Trees.levels(tree))

function levellabel(tree, level)
    nodes = H2Trees.nodesatlevel(tree, level)
    nleaves = count(node -> H2Trees.isleaf(tree, node), nodes)
    return "level $level — $(length(nodes)) nodes, $nleaves leaves"
end

function titleannotation(level)
    return [
        attr(;
            text=levellabel(tree, level),
            showarrow=false,
            xref="paper",
            yref="paper",
            x=0.5,
            y=1.0,
        ),
    ]
end

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
            line=attr(; color=0:(length(x) - 1), colorscale="Viridis", width=4),
            marker=attr(; size=3, color=0:(length(x) - 1), colorscale="Viridis"),
            customdata=ids,
            hovertemplate="node %{customdata}<br>(%{x:.2f}, %{y:.2f}, %{z:.2f})<extra></extra>",
            visible=(i == 1),
        ),
    )
end

steps = [
    attr(;
        label="$level",
        method="update",
        args=[
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
        # Fixed cubic aspect: deeper levels occupy a tiny part of the root box, and without this
        # each level would be rescaled to fill the frame, hiding exactly that.
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
