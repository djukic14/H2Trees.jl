using H2Trees
using PlotlyJS
using StaticArrays

# One point per finest cell, so the tree is complete and the traversal is the full curve.
side = 16
points = [SVector(i + 0.5) for i in 0:(side - 1)]

tree = buildtree(points; builder=TwoNTreeBuilder(; minhalfsize=0.5, minvalues=0))

# `nodesatlevel` is already Hilbert-ordered within each level. In 1D that is left to right.
levelcurves = map(H2Trees.levels(tree)) do level
    nodes = H2Trees.nodesatlevel(tree, level)
    return ([H2Trees.center(tree, node)[1] for node in nodes], nodes)
end

# --- hide-from-docs ---
treelevels = collect(H2Trees.levels(tree))

# Every cell of a level in ONE trace, separated by NaN (see the 2D figure for why): the slider
# toggles trace visibility, so fewer traces means a cheaper mask and smaller emitted HTML.
function cellbars(tree, level)
    xs = Float64[]
    ys = Float64[]
    for node in H2Trees.nodesatlevel(tree, level)
        cx = H2Trees.center(tree, node)[1]
        hs = H2Trees.halfsize(tree, node)
        append!(xs, [cx - hs, cx + hs, NaN])
        append!(ys, [0.0, 0.0, NaN])
    end
    return xs, ys
end

traces = GenericTrace[]
levelof = Int[]
for (i, level) in enumerate(treelevels)
    x, ids = levelcurves[i]
    bx, by = cellbars(tree, level)

    push!(
        traces,
        scatter(;
            x=bx,
            y=by,
            mode="lines",
            line=attr(; color="rgb(225,225,225)", width=10),
            hoverinfo="skip",
            visible=(i == 1),
        ),
    )
    push!(levelof, level)

    push!(
        traces,
        scatter(;
            x=x,
            y=zeros(length(x)),
            mode="lines+markers",
            line=attr(; color="rgb(150,150,150)", width=1.5),
            marker=attr(;
                size=max(5, 11 - i), color=0:(length(x) - 1), colorscale="Viridis"
            ),
            customdata=ids,
            hovertemplate="node %{customdata}<br>center %{x}<extra></extra>",
            visible=(i == 1),
        ),
    )
    push!(levelof, level)
end

function levellabel(level)
    return "level $level — $(length(H2Trees.nodesatlevel(tree, level))) node" *
           (length(H2Trees.nodesatlevel(tree, level)) == 1 ? "" : "s")
end

function titleannotation(level)
    return [
        attr(;
            text=levellabel(level),
            showarrow=false,
            xref="paper",
            yref="paper",
            x=0.5,
            y=1.1,
        ),
    ]
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
        margin=attr(; l=20, r=20, t=50, b=30),
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis=attr(; title="x", zeroline=false),
        # The single row carries no information; level is shown by the slider.
        yaxis=attr(; visible=false, range=[-1, 1]),
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
