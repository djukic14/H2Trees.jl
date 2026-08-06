using CompScienceMeshes
using H2Trees
using PlotlyJS

m = meshsphere(1.0, 0.1)
tree = TwoNTree(vertices(m); builder=TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0))

traces = [wireframe(skeleton(m, 1))]

for node in H2Trees.LevelIterator(tree, 4)
    push!(traces, H2Trees.tracecube(tree, node; mode="lines", line_color=:pink))
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
