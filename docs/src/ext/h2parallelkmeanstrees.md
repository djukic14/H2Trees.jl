# KMeansTree

`KMeansTree(points, numberofclusters; ...)` builds a `BoundingBallTree` by recursively
partitioning points with k-means.

Construction outline:

1. A root cluster is fitted and the root center/radius are initialized.
2. Each node is split with k-means into `numberofclusters` child clusters while
   `length(values(node)) >= max(minvalues, numberofclusters)`.
3. Children store point indices and their local center/radius.
4. Internal-node point lists are emptied after splitting, and node ordering is adjusted.
5. Node radii are finalized via `updateradii!` (default: `boundingsphere`).

!!! warning

    If you use `updateradii=H2Trees.unsafemaxradiusboundingsphere` or
    `updateradii=H2Trees.noboundingsphereupdate`, the radii are not reliably updated.
    This can lead to incorrect results when using iterators in H2Trees.
    Proceed at your own risk.

```@example kmeanstree_1
using CompScienceMeshes
using H2Trees
using PlotlyJS
using ParallelKMeans

m = meshsphere(1.0, 0.1)
tree = KMeansTree(vertices(m), 4; minvalues=60)

traces = [wireframe(skeleton(m, 1))] #hide

for node in H2Trees.LevelIterator(tree, 4) #hide
    push!( #hide
        traces, #hide
        H2Trees.traceball(tree, node; colorscale=[[0, :pink], [1, :pink]],  #hide
        opacity=0.6, showscale=false), #hide
    ) #hide
end #hide

p = PlotlyJS.plot( #hide
    traces, #hide
    Layout(; #hide
        scene=attr(; #hide
            xaxis=attr(; visible=false), #hide
            yaxis=attr(; visible=false), #hide
            zaxis=attr(; visible=false), #hide
        ), #hide
    ), #hide
) #hide

savefig(p, "kmeans_traceball_1.html"); #hide
nothing #hide
```

```@raw html
<object data="../kmeans_traceball_1.html" type="text/html"  style="width:100%; height:50vh;"> </object>
```

```@example kmeanstree_2
using Logging #hide
using CompScienceMeshes
using H2Trees
using PlotlyJS
using ParallelKMeans

m = meshsphere(1.0, 0.1)
with_logger(NullLogger()) do #hide
tree = KMeansTree(
    vertices(m),
    4;
    minvalues=60,
    updateradii=H2Trees.noboundingsphereupdate,
)

traces = [wireframe(skeleton(m, 1))] #hide

for node in H2Trees.LevelIterator(tree, 4) #hide
    push!( #hide
        traces, #hide
        H2Trees.traceball(tree, node; colorscale=[[0, :pink], [1, :pink]],  #hide
        opacity=0.6, showscale=false), #hide
    ) #hide
end #hide

p = PlotlyJS.plot( #hide
    traces, #hide
    Layout(; #hide
        scene=attr(; #hide
            xaxis=attr(; visible=false), #hide
            yaxis=attr(; visible=false), #hide
            zaxis=attr(; visible=false), #hide
        ), #hide
    ), #hide
) #hide

savefig(p, "kmeans_traceball_2.html"); #hide
end #hide
nothing #hide
```

```@raw html
<object data="../kmeans_traceball_2.html" type="text/html"  style="width:100%; height:50vh;"> </object>
```

```@example kmeanstree_3
using Logging #hide
using CompScienceMeshes
using H2Trees
using PlotlyJS
using ParallelKMeans

m = meshsphere(1.0, 0.1)
with_logger(NullLogger()) do #hide
tree = KMeansTree(
    vertices(m),
    4;
    minvalues=60,
    updateradii=H2Trees.unsafemaxradiusboundingsphere
)

traces = [wireframe(skeleton(m, 1))] #hide

for node in H2Trees.LevelIterator(tree, 4) #hide
    push!( #hide
        traces, #hide
        H2Trees.traceball(tree, node; colorscale=[[0, :pink], [1, :pink]],  #hide
        opacity=0.6, showscale=false), #hide
    ) #hide
end #hide

p = PlotlyJS.plot( #hide
    traces, #hide
    Layout(; #hide
        scene=attr(; #hide
            xaxis=attr(; visible=false), #hide
            yaxis=attr(; visible=false), #hide
            zaxis=attr(; visible=false), #hide
        ), #hide
    ), #hide
) #hide

savefig(p, "kmeans_traceball_3.html"); #hide
end #hide
nothing #hide
```

```@raw html
<object data="../kmeans_traceball_3.html" type="text/html"  style="width:100%; height:50vh;"> </object>
```
