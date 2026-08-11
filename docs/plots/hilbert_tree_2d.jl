using H2Trees
using PlotlyJS
using StaticArrays

# A complete 8x8 grid, one point per finest cell, so every box on every level is occupied and
# the traversal below is the full Hilbert curve rather than a subsequence of it.
side = 8
points = [SVector(i + 0.5, j + 0.5) for i in 0:(side - 1) for j in 0:(side - 1)]

tree = buildtree(points; builder=TwoNTreeBuilder(; minhalfsize=0.5, minvalues=0))

# `nodesatlevel` is Hilbert-ordered within each level, so node centers trace the curve directly.
levelcurves = map(H2Trees.levels(tree)) do level
    nodes = H2Trees.nodesatlevel(tree, level)
    centers = [H2Trees.center(tree, node) for node in nodes]
    return (getindex.(centers, 1), getindex.(centers, 2), nodes)
end

# --- hide-from-docs ---
treelevels = collect(H2Trees.levels(tree))

# Every box of a level in ONE trace, separated by NaN, rather than one trace per box: the
# slider below toggles trace visibility, and a mask over 2 traces per level is far cheaper (and
# far smaller in the emitted HTML) than one over ~64.
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
            marker=attr(;
                size=max(4, 9 - i), color=0:(length(x) - 1), colorscale="Viridis"
            ),
            customdata=ids,
            hovertemplate="node %{customdata}<br>(%{x}, %{y})<extra></extra>",
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
            attr(;
                annotations=[
                    attr(;
                        text="level $level — $(length(H2Trees.nodesatlevel(tree, level))) nodes",
                        showarrow=false,
                        xref="paper",
                        yref="paper",
                        x=0.5,
                        y=1.05,
                    ),
                ],
            ),
        ],
    ) for level in treelevels
]

p = PlotlyJS.plot(
    traces,
    Layout(;
        showlegend=false,
        margin=attr(; l=20, r=20, t=40, b=20),
        plot_bgcolor="rgba(0,0,0,0)",
        # Square cells: without `scaleanchor` the curve is sheared and the U-shapes are hard to
        # recognise.
        xaxis=attr(; visible=false, scaleanchor="y", scaleratio=1),
        yaxis=attr(; visible=false),
        annotations=[
            attr(;
                text="level $(first(treelevels)) — $(length(H2Trees.nodesatlevel(tree, first(treelevels)))) nodes",
                showarrow=false,
                xref="paper",
                yref="paper",
                x=0.5,
                y=1.05,
            ),
        ],
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
