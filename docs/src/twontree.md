# TwoNTree

The [`TwoNTree`](@ref) is a $2^n$-tree for organizing points in $\mathbb{R}^n$.
For 3D points, this is an octree.

You can build a [`TwoNTree`](@ref) from a set of points and a minimum `halfsize`.
The most important options are:

- `minvalues`: minimum number of points required for a box before it is subdivided.
- `maxprotrusion`: maximum allowed protrusion of an element outside its box, normalized by `2*halfsize`.
- `computeprotrusion`: functor used to compute protrusion.

For example, `maxprotrusion=0.5` means no element may protrude by more than `halfsize`.

By default, `computeprotrusion=ComputeProtrusionFunctor()`, which treats elements as points and therefore gives zero protrusion.
For more complex elements, provide a custom protrusion functor.

In the Galerkin case, the tree can be constructed as follows

```@example TwoNTree1
using  CompScienceMeshes # hide
using H2Trees # hide

m = meshsphere(1.0, 0.1)
tree = TwoNTree(vertices(m), 0.1; minvalues=60, maxprotrusion=0.5)
```

Alternatively, the tree can be constructed with a minimum `halfsize` of 0:

```@example TwoNTree2
using  CompScienceMeshes # hide
using H2Trees # hide

m = meshsphere(1.0, 0.1)
tree = TwoNTree(vertices(m), 0.0; minvalues=60)
```

In the Petrov-Galerkin case, the tree can be constructed by providing two sets of points

```@example TwoNTree3
using CompScienceMeshes
using H2Trees

mx = meshsphere(1.0, 0.1)
my = meshsphere(2.0, 0.1)

tree = TwoNTree(vertices(mx), vertices(my), 0.1)
```

This creates a [`BlockTree`](@ref) with two [`TwoNTree`](@ref)s.
!!! note
    The trees are configured such that both trees have the same `halfsize` at the same `level`.
    This means that not every tree starts at `level` 1.
