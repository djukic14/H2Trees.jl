# Visualizing a TwoNTree

```@example plot_twontree
using CompScienceMeshes
using H2Trees
using PlotlyJS

m = meshsphere(1.0, 0.1)
tree = TwoNTree(vertices(m), 0.1)

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

savefig(p, "sphere_tracecube.html"); # hide
nothing #hide
```

```@raw html
<object data="../sphere_tracecube.html" type="text/html"  style="width:100%; height:50vh;"> </object>
```
