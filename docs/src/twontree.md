# TwoNTree

The [`TwoNTree`](@ref) is a 2ⁿ-tree for organizing points in ℝⁿ.
In the case of three-dimensional points, this results in an octree.
The [`TwoNTree`](@ref) can be constructed from a set of points and a minimum `halfsize`.
Additionally, the `minvalues` parameter can be set to control the subdivision of the tree.
A box is only further subdivided if it contains at least `minvalues` points.
The parameter `maxprotrusion` can be set such that no element in the tree protrudes more than `maxprotrusion` out of the box containing it.
The protrusion is normalized with regard to the 2*`halfsize` of the box, so that `maxprotrusion=0.5` means that no element protrudes more than half of the 2*`halfsize` out of the box.
The computation of the protrusion defaults to `computeprotrusion=ComputeProtrusionFunctor()` which treats the elements of the tree as points, so that the protrusion is zero for all elements.
For more complex elements, a custom protrusion functor can be provided that computes the protrusion of the element with respect to the box.

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
