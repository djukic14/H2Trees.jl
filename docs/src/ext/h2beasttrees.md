# BEAST

The basis functions defined in [BEAST.jl](https://github.com/krcools/BEAST.jl) can be sorted into trees.

## TwoNTree together with BEAST

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
