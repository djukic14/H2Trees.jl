# Forest

A [`Forest`](@ref) stores multiple trees and provides a simple, array-like interface
for accessing and iterating over them.

In practice, forests are useful when your geometry is disconnected (like in a multibody problem):
each connected component becomes its own tree, and all trees are grouped into one
container.

## Using a Forest

The example below focuses on how to use a forest once it is available.

```@example Forest1
using Metis, CompScienceMeshes, BEAST # hide
using H2Trees # hide

m = CompScienceMeshes.readmesh( # hide
    joinpath(pkgdir(H2Trees), "test", "assets", "in", "twospheres.in") # hide
) # hide

X = lagrangecxd0(m) # hide

forest = H2Trees.MetisForest(X, 4; updateradii=H2Trees.boundingsphere) # hide

@show length(forest)

for tree in forest
    println("Tree with ", H2Trees.numberofnodes(tree), " nodes")
end
```

`length(forest)` returns the number of component trees, and iterating over
`forest` yields each tree in sequence.

## Working with `Forest`

`Forest` supports the usual collection operations:

- `length(forest)` gives the number of trees.
- `forest[i]` accesses the `i`-th tree.
- `for tree in forest` iterates through all trees.

More forest-related functionality and convenience utilities will be added in future releases.
