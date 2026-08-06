using CompScienceMeshes, PlotlyJS
using H2Trees
using ParallelKMeans

m = meshsphere(1.0, 0.1)
boxtree = TwoNTree(vertices(m); builder=TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0))
balltree = KMeansTree(
    vertices(m); builder=KMeansTreeBuilder(; numberofclusters=4, minvalues=60)
)

# --- hide-from-docs ---
traces = [wireframe(skeleton(m, 1))]
for node in H2Trees.LevelIterator(boxtree, 4)
    push!(traces, H2Trees.tracecube(boxtree, node; mode="lines", line_color=:pink))
end
for node in H2Trees.LevelIterator(balltree, 3)
    push!(
        traces,
        H2Trees.traceball(
            balltree,
            node;
            colorscale=[[0, :blue], [1, :blue]],
            opacity=0.25,
            showscale=false,
        ),
    )
end

p = PlotlyJS.plot(
    traces,
    Layout(;
        scene=attr(;
            xaxis=attr(; visible=false),
            yaxis=attr(; visible=false),
            zaxis=attr(; visible=false),
        ),
        showlegend=false,
    ),
)
