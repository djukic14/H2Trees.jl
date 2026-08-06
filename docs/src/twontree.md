# TwoNTree

The [`TwoNTree`](@ref) is a $2^n$-tree for organizing points in $\mathbb{R}^n$.
For 3D points, this is an octree.

You can build a [`TwoNTree`](@ref) from a set of points through the builder-based construction API.
The most important [`TwoNTreeBuilder`](@ref) options are:

- `minhalfsize`: minimum half-size of a leaf box. The default `0` disables this stopping criterion.
- `minvalues`: minimum number of points required for a box before it is subdivided. The default is `70`.
- `protrusion`: protrusion policy, see [Protrusion Policy](protrusion.md).

In the Galerkin case, the tree can be constructed as follows

```@example TwoNTree1
using  CompScienceMeshes # hide
using H2Trees # hide

m = meshsphere(1.0, 0.1)
tree = H2Trees.buildtree(
    vertices(m);
    builder=TwoNTreeBuilder(;
        minhalfsize=0.1,
        minvalues=60,
        protrusion=H2Trees.ProtrusionCheck(; max=0.5),
    ),
)
```

Alternatively, the tree can be constructed with a minimum `halfsize` of 0:

```@example TwoNTree2
using  CompScienceMeshes # hide
using H2Trees # hide

m = meshsphere(1.0, 0.1)
tree = H2Trees.buildtree(vertices(m); builder=TwoNTreeBuilder(; minhalfsize=0.0, minvalues=60))
```

In the Petrov-Galerkin case, two sets of points build a [`BlockTree`](@ref) of two
[`TwoNTree`](@ref)s instead — see [BlockTree and Petrov Trees](blocktree_petrov.md) for
construction and the level-scale invariant both sides must satisfy.

Internally, `TwoNTree` stores cached topology in a tree index; see
[Tree Access and Values](tree_access.md) for the accessor API.
