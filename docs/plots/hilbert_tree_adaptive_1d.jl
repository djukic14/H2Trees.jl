using H2Trees
using PlotlyJS
using StaticArrays
using Random

# Dense cluster plus a sparse tail: `minvalues` keeps splitting the cluster long after the tail
# has stopped, so the tree runs far deeper than the uniform figure above and its leaves sit on
# many different levels.
Random.seed!(21)
points = vcat([SVector(0.02 * randn()) for _ in 1:300], [SVector(randn()) for _ in 1:60])

tree = buildtree(points; builder=TwoNTreeBuilder(; minvalues=8))

levelcurves = map(H2Trees.levels(tree)) do level
    nodes = H2Trees.nodesatlevel(tree, level)
    return ([H2Trees.center(tree, node)[1] for node in nodes], nodes)
end

# --- hide-from-docs ---
treelevels = collect(H2Trees.levels(tree))

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
            y=1.1,
        ),
    ]
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
            marker=attr(; size=7, color=0:(length(x) - 1), colorscale="Viridis"),
            customdata=ids,
            hovertemplate="node %{customdata}<br>center %{x:.4f}<extra></extra>",
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
        margin=attr(; l=20, r=20, t=50, b=30),
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis=attr(; title="x", zeroline=false),
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
