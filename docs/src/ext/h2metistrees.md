# Metis

The METIS extension provides graph-based partitioning for tree construction.
Once `Metis.jl` is loaded, H2Trees can use `metispartition` internally and expose
high-level constructors such as [`MetisTree`](@ref) and [`MetisForest`](@ref).

## Adjacency Graph and Weights

For BEAST spaces, the graph and weights are typically built with
`g, w = H2Trees.adjacencygraph(X)`.

- `g` is a graph on basis functions.
- A graph edge between two basis functions means their supports touch through the
  underlying mesh.
- `w` contains one vertex weight per basis function and is assembled from mesh
  element areas distributed over all basis functions supported on each element.

These weights are then used by METIS to bias partitioning.

```@example metis_graph
using H2Trees
using Metis, Graphs
using BEAST, CompScienceMeshes

m = meshicosphere(1)
X = lagrangecxd0(m)

g, w = H2Trees.adjacencygraph(X)
```

Set `splitunconnectedpartitions=true` when you want each output partition to be
post-processed into connected components.

## Visualizing a `MetisTree`

A tree node can be visualized by coloring basis functions that belong to each child cluster.

```@example metis_tree_vis
using H2Trees
using Metis
using BEAST
using CompScienceMeshes
using PlotlyJS

m = meshicosphere(10)
X = lagrangecxd0(m)

tree = H2Trees.MetisTree(X, 4)

childcolors = [
    "rgb(0.5490196078431373, 0.00784313725490196, 0.45098039215686275)",
    "rgb(0.6, 0.38823529411764707, 0.18823529411764706)",
    "rgb(0.5019607843137255, 0.7725490196078432, 0.37254901960784315)",
    "rgb(0.7019607843137254, 0.9490196078431372, 0.9921568627450981)",
]

color = ["rgba(0,0,0,0)" for _ in 1:numcells(m)]
traces = [wireframe(skeleton(m, 1); width=1)] # hide
push!(traces, mesh3d(m; facecolor=color)) # hide

node = H2Trees.root(tree)
for (i, child) in enumerate(H2Trees.children(tree, node))
    for val in H2Trees.values(tree, child)
        for fn in X.fns[val]
            color[fn.cellid] = childcolors[i]
        end
    end
end

p = PlotlyJS.plot( # hide
    traces, # hide
    Layout(; # hide
        scene=attr(; # hide
            aspectmode="data", # hide
                xaxis=attr(; visible=false), # hide
                yaxis=attr(; visible=false), # hide
                zaxis=attr(; visible=false), # hide
        ), # hide
    ), # hide
) # hide

savefig(p, "metis_tree_clusters.html"); # hide
nothing # hide
```

```@raw html
<object data="../metis_tree_clusters.html" type="text/html"  style="width:100%; height:50vh;"> </object>
```

## `MetisForest` for Disconnected Geometry

`MetisForest` is useful when your adjacency graph has multiple connected
components. One tree is created per connected component and wrapped in a
[`Forest`](@ref).

```@example metis_forest_vis
using H2Trees
using Metis
using BEAST
using CompScienceMeshes
using PlotlyJS

m = CompScienceMeshes.readmesh(
    joinpath(pkgdir(H2Trees), "test", "assets", "in", "twospheres2.in")
)
X = lagrangecxd0(m)

forest = H2Trees.MetisForest(X, 4)

childcolors = [
    "rgb(0.5490196078431373, 0.00784313725490196, 0.45098039215686275)",
    "rgb(0.6, 0.38823529411764707, 0.18823529411764706)",
    "rgb(0.5019607843137255, 0.7725490196078432, 0.37254901960784315)",
    "rgb(0.7019607843137254, 0.9490196078431372, 0.9921568627450981)",
    "rgb(0.5725490196078431, 0.16470588235294117, 0.34901960784313724)",
    "rgb(0.615686274509804, 0.5137254901960784, 0.11764705882352941)",
    "rgb(0.4, 0.8470588235294118, 0.611764705882353)",
    "rgb(0.4235294117647059, 0.9215686274509803, 0.8588235294117647)",
]

color = ["rgba(0,0,0,0)" for _ in 1:numcells(m)]
traces = [wireframe(skeleton(m, 1); width=1)] # hide
push!(traces, mesh3d(m; facecolor=color)) # hide

node = H2Trees.root(forest[1])
for (i, child) in enumerate(H2Trees.children(forest[1], node))
    for val in H2Trees.values(forest[1], child)
        for fn in X.fns[val]
            color[fn.cellid] = childcolors[i]
        end
    end
end

node = H2Trees.root(forest[2])
for (i, child) in enumerate(H2Trees.children(forest[2], node))
    for val in H2Trees.values(forest[2], child)
        for fn in X.fns[val]
            color[fn.cellid] = childcolors[i + 4]
        end
    end
end

p = PlotlyJS.plot( # hide
    traces, # hide
    Layout(; # hide
        scene=attr(; # hide
        aspectmode="data", # hide
        xaxis=attr(; visible=false), # hide
        yaxis=attr(; visible=false), # hide
        zaxis=attr(; visible=false), # hide
        ), # hide
    ), # hide
) # hide

savefig(p, "metis_forest_clusters.html"); # hide
nothing # hide
```

```@raw html
<object data="../metis_forest_clusters.html" type="text/html"  style="width:100%; height:50vh;"> </object>
```
