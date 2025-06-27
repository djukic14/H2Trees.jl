# BEAST

The basis functions defined in [BEAST.jl](https://github.com/krcools/BEAST.jl) can be sorted into trees.

## TwoNTree

For [`TwoNTree`](@ref)s this can, for example, be done like this for the Galerkin case

```@example BEAST1
using BEAST, CompScienceMeshes
using H2Trees

m = meshsphere(1.0, 0.1)
X = raviartthomas(m)
tree = TwoNTree(X, 0.1)
```

and for the Petrov-Galerkin case

```@example BEAST2
using BEAST, CompScienceMeshes
using H2Trees
using PlotlyJS

m = meshsphere(1.0, 0.1)
X = raviartthomas(m)
Y = buffachristiansen(m)
tree = TwoNTree(X, Y, 0.1)
```

## QuadPointsTree

Similarly [`QuadPointsTree`](@ref)s can be constructed for [BEAST.jl](https://github.com/krcools/BEAST.jl) spaces for the Galerkin case

```@example BEAST3
using BEAST, CompScienceMeshes
using H2Trees

m = meshsphere(1.0, 0.1)
X = raviartthomas(m)
tree = QuadPointsTree(X, 0.1)
```

and for the Petrov-Galerkin case

```@example BEAST4
using BEAST, CompScienceMeshes
using H2Trees
using PlotlyJS

m = meshsphere(1.0, 0.1)
X = raviartthomas(m)
Y = buffachristiansen(m)
tree = QuadPointsTree(X, Y, 0.1)
```
