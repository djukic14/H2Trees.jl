# Metis

The METIS extension provides graph-based partitioning for tree construction.
Once `Metis.jl` is loaded, H2Trees can use `metispartition` internally and expose
high-level constructors such as [`MetisTree`](@ref) and [`MetisForest`](@ref).

Pass explicit builders. [`MetisTreeBuilder`](@ref) and [`MetisForestBuilder`](@ref) store the
partitioning options explicitly.

```julia
tree = MetisTree(
    X;
    builder=MetisTreeBuilder(;
        numdivisions=4,
    ),
)

forest = MetisForest(
    X;
    builder=MetisForestBuilder(;
        treebuilder=MetisTreeBuilder(;
            numdivisions=4,
        ),
    ),
)
```

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

```@eval
using H2Trees
include(joinpath(pkgdir(H2Trees), "docs", "plotutils.jl"))
displayedcode(joinpath(pkgdir(H2Trees), "docs", "plots", "metis_tree_clusters.jl"))
```

```@raw html
<object data="../../assets/plots/metis_tree_clusters.html" type="text/html"  style="width:100%; height:50vh;"> </object>
```

## `MetisForest` for Disconnected Geometry

`MetisForest` is useful when your adjacency graph has multiple connected
components. One tree is created per connected component and wrapped in a
[`Forest`](@ref).

```@eval
using H2Trees
include(joinpath(pkgdir(H2Trees), "docs", "plotutils.jl"))
displayedcode(joinpath(pkgdir(H2Trees), "docs", "plots", "metis_forest_clusters.jl"))
```

```@raw html
<object data="../../assets/plots/metis_forest_clusters.html" type="text/html"  style="width:100%; height:50vh;"> </object>
```
