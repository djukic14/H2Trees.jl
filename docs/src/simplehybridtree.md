# Simple Hybrid Tree

A [`SimpleHybridTree`](@ref) is a tree that is split into two parts: an upper tree and a lower tree.
The upper tree spans the levels from `minlevel` to `hybridlevel`, while the lower tree spans the levels from `hybridlevel`+1 to `maxlevel`.
This allows the 𝓗² Method to treat both parts of the tree differently.

A [`SimpleHybridTree`](@ref) can be constructed from a [`TwoNTree`](@ref) by specifying the `hybridhalfsize` parameter.
The `hybridhalfsize` determines the level at which the tree is split.

In the Galerkin case, a hybrid [`TwoNTree`](@ref) can be constructed as follows:

```@example TwoNTree1
using  CompScienceMeshes # hide
using H2Trees # hide

m = meshsphere(1.0, 0.1)
tree = SimpleHybridTree(TwoNTree(vertices(m), 0.1); hybridhalfsize=0.2)
```

Translating plans can be split into two with the `splitplan` function.
