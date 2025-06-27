# PlotlyJS

[PlotlyJS.jl](https://github.com/JuliaPlots/PlotlyJS.jl) can be used to visualize the clusters of a tree.

## Visualizing a TwoNTree

For [`TwoNTree`](@ref)s we have the helper-function [`tracecube`](@ref), which can, for example, used like this

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

## Visualizing a BoundingBallTree

For [`BoundingBallTree`](@ref)s we have the helper-function [`traceball`](@ref), which can, for example, used like this

```@example plot_boundingball
using CompScienceMeshes
using H2Trees
using PlotlyJS
using ParallelKMeans

m = meshsphere(1.0, 0.1)
tree = KMeansTree(vertices(m), 4; minvalues=60)

traces = [wireframe(skeleton(m, 1))]

for node in H2Trees.LevelIterator(tree, 4)
    push!(
        traces,
        H2Trees.traceball(tree, node; colorscale=[[0, :pink], [1, :pink]], 
        opacity=0.6, showscale=false),
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
    ),
)

savefig(p, "sphere_traceball.html"); # hide
nothing #hide
```

```@raw html
<object data="../sphere_traceball.html" type="text/html"  style="width:100%; height:50vh;"> </object>
```
