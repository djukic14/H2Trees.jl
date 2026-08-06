# PlotlyJS

[PlotlyJS.jl](https://github.com/JuliaPlots/PlotlyJS.jl) can be used to visualize the clusters of a tree.

## Visualizing a TwoNTree

For [`TwoNTree`](@ref)s we have the helper-function [`tracecube`](@ref), which can, for example, used like this

```@eval
using H2Trees
include(joinpath(pkgdir(H2Trees), "docs", "plotutils.jl"))
displayedcode(joinpath(pkgdir(H2Trees), "docs", "plots", "sphere_tracecube.jl"))
```

```@raw html
<object data="../../assets/plots/sphere_tracecube.html" type="text/html"  style="width:100%; height:50vh;"> </object>
```

## Visualizing a BoundingBallTree

For [`BoundingBallTree`](@ref)s we have the helper-function [`traceball`](@ref), which can, for example, used like this

```@eval
using H2Trees
include(joinpath(pkgdir(H2Trees), "docs", "plotutils.jl"))
displayedcode(joinpath(pkgdir(H2Trees), "docs", "plots", "sphere_traceball.jl"))
```

```@raw html
<object data="../../assets/plots/sphere_traceball.html" type="text/html"  style="width:100%; height:50vh;"> </object>
```
