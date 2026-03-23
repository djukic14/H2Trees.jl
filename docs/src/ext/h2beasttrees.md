# BEAST

The basis functions defined in [BEAST.jl](https://github.com/krcools/BEAST.jl) can be sorted into trees.

## TwoNTree together with BEAST

For basis functoins defined in [BEAST.jl](https://github.com/krcools/BEAST.jl) the computation of the protrusion defaults to `BEASTProtrusionFunctor` which means that setting the `maxprotrusion` is actually meaningfull.
For [`TwoNTree`](@ref)s this can, for example, be done like this for the Galerkin case

```@example BEAST1
using BEAST, CompScienceMeshes
using H2Trees

m = meshsphere(1.0, 0.1)
X = raviartthomas(m)
tree = TwoNTree(X, 0.0; maxprotrusion=0.3, minvalues=200)
```

and for the Petrov-Galerkin case

```@example BEAST2
using BEAST, CompScienceMeshes
using H2Trees
using PlotlyJS

m = meshsphere(1.0, 0.1)
X = raviartthomas(m)
Y = buffachristiansen(m)
tree = TwoNTree(X, Y, 0.0; testmaxprotrusion=0.3, trialmaxprotrusion=0.3)
```
