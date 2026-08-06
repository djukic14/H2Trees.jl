# KMeansTree

`KMeansTree(points; builder=KMeansTreeBuilder(...))` builds a `BoundingBallTree` by recursively
partitioning points with k-means.

Construction outline:

1. A root cluster is fitted and the root center/radius are initialized.
2. Each node is split with k-means into `numberofclusters` child clusters while
   `length(values(node)) >= max(minvalues, numberofclusters)`.
3. Children store point indices and their local center/radius.
4. Internal-node point lists are emptied after splitting, and node ordering is adjusted.
5. Node radii are finalized via `updateradii!` (default: `boundingsphere`).

Pass a [`KMeansTreeBuilder`](@ref), which stores all construction settings explicitly:

```julia
tree = KMeansTree(
    points;
    builder=KMeansTreeBuilder(;
        numberofclusters=numberofclusters,
        minvalues=60,
    ),
)
```

!!! warning

    If you use `updateradii=H2Trees.unsafemaxradiusboundingsphere` or
    `updateradii=H2Trees.noboundingsphereupdate`, the radii are not reliably updated.
    This can lead to incorrect results when using iterators in H2Trees.
    Proceed at your own risk.

```@eval
using H2Trees
include(joinpath(pkgdir(H2Trees), "docs", "plotutils.jl"))
displayedcode(joinpath(pkgdir(H2Trees), "docs", "plots", "kmeans_traceball_1.jl"))
```

```@raw html
<object data="../../assets/plots/kmeans_traceball_1.html" type="text/html"  style="width:100%; height:50vh;"> </object>
```

```@eval
using H2Trees
include(joinpath(pkgdir(H2Trees), "docs", "plotutils.jl"))
displayedcode(joinpath(pkgdir(H2Trees), "docs", "plots", "kmeans_traceball_2.jl"))
```

```@raw html
<object data="../../assets/plots/kmeans_traceball_2.html" type="text/html"  style="width:100%; height:50vh;"> </object>
```

```@eval
using H2Trees
include(joinpath(pkgdir(H2Trees), "docs", "plotutils.jl"))
displayedcode(joinpath(pkgdir(H2Trees), "docs", "plots", "kmeans_traceball_3.jl"))
```

```@raw html
<object data="../../assets/plots/kmeans_traceball_3.html" type="text/html"  style="width:100%; height:50vh;"> </object>
```
