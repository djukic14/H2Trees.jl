using H2Trees
using PlotlyJS
using StaticArrays
using Random

# Strongly non-uniform: a dense core plus a sparse halo. With `minvalues` stopping subdivision
# once a box holds few enough points, the core keeps splitting long after the halo has stopped,
# so the tree is deep, adaptive, and has leaves on many different levels.
Random.seed!(22)
points = vcat(
    [SVector((0.03 .* randn(2))...) for _ in 1:600], [SVector(randn(2)...) for _ in 1:150]
)

tree = buildtree(points; builder=TwoNTreeBuilder(; minvalues=12))

# Same walk as the uniform figure: `nodesatlevel` is increasing node id, which within a level is
# Hilbert order. Only the geometry changed.
levelcurves = map(H2Trees.levels(tree)) do level
    nodes = H2Trees.nodesatlevel(tree, level)
    centers = [H2Trees.center(tree, node) for node in nodes]
    return (getindex.(centers, 1), getindex.(centers, 2), nodes)
end

# --- hide-from-docs ---
treelevels = collect(H2Trees.levels(tree))

function boxoutlines(tree, level)
    xs = Float64[]
    ys = Float64[]
    for node in H2Trees.nodesatlevel(tree, level)
        cx, cy = H2Trees.center(tree, node)
        hs = H2Trees.halfsize(tree, node)
        append!(xs, [cx - hs, cx + hs, cx + hs, cx - hs, cx - hs, NaN])
        append!(ys, [cy - hs, cy - hs, cy + hs, cy + hs, cy - hs, NaN])
    end
    return xs, ys
end

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
            y=1.05,
        ),
    ]
end

traces = GenericTrace[]
levelof = Int[]
for (i, level) in enumerate(treelevels)
    x, y, ids = levelcurves[i]
    bx, by = boxoutlines(tree, level)

    push!(
        traces,
        scatter(;
            x=bx,
            y=by,
            mode="lines",
            line=attr(; color="rgb(225,225,225)", width=1),
            hoverinfo="skip",
            visible=(i == 1),
        ),
    )
    push!(levelof, level)

    push!(
        traces,
        scatter(;
            x=x,
            y=y,
            mode="lines+markers",
            line=attr(; color="rgb(150,150,150)", width=1.5),
            marker=attr(; size=6, color=0:(length(x) - 1), colorscale="Viridis"),
            customdata=ids,
            hovertemplate="node %{customdata}<br>(%{x:.3f}, %{y:.3f})<extra></extra>",
            visible=(i == 1),
        ),
    )
    push!(levelof, level)
end

steps = [
    attr(;
        label="$level",
        method="update",
        args=[
            attr(; visible=[l == level for l in levelof]),
            attr(; annotations=titleannotation(level)),
        ],
    ) for level in treelevels
]

p = PlotlyJS.plot(
    traces,
    Layout(;
        showlegend=false,
        margin=attr(; l=20, r=20, t=40, b=20),
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis=attr(; visible=false, scaleanchor="y", scaleratio=1),
        yaxis=attr(; visible=false),
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
