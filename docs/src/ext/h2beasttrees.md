# BEAST

The basis functions defined in [BEAST.jl](https://github.com/krcools/BEAST.jl) can be sorted into trees.

## TwoNTree together with BEAST

For BEAST spaces, `TwoNTree(space)` and `buildtree(space)` use the builder workflow and resolve the default `AutoProtrusionCheck()` to `ProtrusionCheck(; max=0.25, compute=BEASTProtrusionFunctor(space))`.
Use `NoProtrusionCheck()` to opt out, or pass an explicit `ProtrusionCheck` to choose another threshold.

For the Galerkin case:

```@example BEAST1
using BEAST, CompScienceMeshes
using H2Trees

m = meshsphere(1.0, 0.1)
X = raviartthomas(m)
tree = TwoNTree(X; builder=TwoNTreeBuilder(; minhalfsize=0.0, minvalues=200))
```

For the Petrov-Galerkin case:

```@example BEAST2
using BEAST, CompScienceMeshes
using H2Trees
using PlotlyJS

m = meshsphere(1.0, 0.1)
X = raviartthomas(m)
Y = buffachristiansen(m)
tree = buildtree(
    X,
    Y;
    builder=BlockTreeBuilder(;
        test=TwoNTreeBuilder(; minhalfsize=0.0, protrusion=ProtrusionCheck(; max=0.3)),
        trial=TwoNTreeBuilder(; minhalfsize=0.0, protrusion=ProtrusionCheck(; max=0.3)),
    ),
)
```

The same BEAST-space entry point also accepts the clustering builders:

```@example BEAST3
using BEAST, CompScienceMeshes
using H2Trees
using ParallelKMeans

m = meshsphere(1.0, 0.1)
X = lagrangecxd0(m)
tree = buildtree(X; builder=KMeansTreeBuilder())
```

For METIS-based trees the graph is derived from the space by default, or can be passed explicitly
with `graphweights=(graph, weights)`:

```@example BEAST4
using BEAST, CompScienceMeshes
using H2Trees
using Metis

m = meshsphere(1.0, 0.1)
X = lagrangecxd0(m)
tree = buildtree(X; builder=MetisTreeBuilder())
forest = buildforest(X; builder=MetisForestBuilder())
```
