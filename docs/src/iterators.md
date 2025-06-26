# Iterators

Iterators allow to traverse the tree.

## Depth First

The [`DepthFirstIterator`](@ref) traverses the tree in a depth first manner. If no node is specified the tree is
traversed from the root node.

```@example depthfirst
using CompScienceMeshes # hide
using H2Trees # hide

m = meshsphere(1.0, 0.1)
tree = TwoNTree(vertices(m), 0.1)
println(collect(H2Trees.DepthFirstIterator(tree)))
```

## Parents Upward

The [`ParentUpwardsIterator`](@ref) is an iterator that iterates over all parent nodes of a given node in
a tree until the root is reached. The last node is the node 0.

```@example parentupwards
using CompScienceMeshes # hide
using H2Trees # hide

m = meshsphere(1.0, 0.1)
tree = TwoNTree(vertices(m), 0.1)
println(collect(H2Trees.ParentUpwardsIterator(tree, 373)))
```

## Children

The [`ChildIterator`](@ref) is an iterator over the children of a node in a tree.

```@example children
using CompScienceMeshes # hide
using H2Trees # hide

m = meshsphere(1.0, 0.1)
tree = TwoNTree(vertices(m), 0.1)
println(collect(H2Trees.children(tree, H2Trees.root(tree))))
```

## Leaves

[`leaves`](@ref) returns an iterator over the leaf nodes in the tree, starting from the specified `node`.
If no node is specified the tree is traversed from the root node.

```@example leaves
using CompScienceMeshes # hide
using H2Trees # hide

m = meshsphere(1.0, 0.1)
tree = TwoNTree(vertices(m), 0.1)
println(collect(H2Trees.leaves(tree)))
```

## Level Iterator

[`LevelIterator`](@ref) return an iterator over the nodes at the specified `level` in the `tree`.

```@example leveliterator
using CompScienceMeshes # hide
using H2Trees # hide

m = meshsphere(1.0, 0.1)
tree = TwoNTree(vertices(m), 0.1)
println(collect(H2Trees.LevelIterator(tree,2)))
```

## Same Level Nodes

[`SameLevelIterator`] returns an iterator over the nodes at the same level as `node` in the `tree`.

```@example sameleveliterator
using CompScienceMeshes # hide
using H2Trees # hide

m = meshsphere(1.0, 0.1)
tree = TwoNTree(vertices(m), 0.1)
println(collect(H2Trees.SameLevelIterator(tree,3)))
```

## Near- and Far Nodes

The [`NearNodeIterator`](@ref) returns, in the Galerkin case, an iterator over the nodes in the tree that are at the same level as the specified
`node` and are near to `node`. Two nodes are near if the function `isnear(tree, nodea, nodeb)` evaluates to true.

```@example galerkinnearnodes
using CompScienceMeshes # hide
using H2Trees # hide

m = meshsphere(1.0, 0.1)
tree = TwoNTree(vertices(m), 0.1)
println("node 4 is at level ", H2Trees.level(tree, 4))
println("nodes near to node 4: ", collect(H2Trees.NearNodeIterator(tree, 4)))
```

In the Petrov-Galerkin case the [`NearNodeIterator`](@ref) returns an iterator over the nodes in the `testtree` that are at the same level as the specified `node` in the `trialtree` and are near to `node`.
Two nodes are near if the function `isnear(testtree, trialtree, testnode, trialnode)` evaluates to true.

```@example galerkinnearnodes
using CompScienceMeshes # hide
using H2Trees # hide

mx = meshsphere(1.0, 0.1)
my = meshsphere(2.0, 0.1)

testtree = TwoNTree(vertices(mx), 0.1)
trialtree = TwoNTree(vertices(my), 0.1)
println("trial node 4 is at level ", H2Trees.level(trialtree, 4))
println("nodes near to trial node 4: ", collect(H2Trees.NearNodeIterator(testtree, trialtree, 4)))
```

The [`FarNodeIterator`](@ref) is defined accordingly.

## Well Separated Nodes

## Filtering Nodes

[`NodeFilterIterator`](@ref) returns an iterator over nodes at the same level as `trialnode` in `trialtree`,
for which the function `filter(testtree, trialtree, testnode, trialnode)` return true.
In the case of a tree with the trait `isBlockTree`, it is assumed that `trialnode` is
belonging to the `trialtree`.
For the Galerkin case, we have

```@example filteringnodes1
using CompScienceMeshes # hide
using H2Trees # hide

m = meshsphere(1.0, 0.1)
tree = TwoNTree(vertices(m), 0.1)
println(collect(H2Trees.NodeFilterIterator(tree, 3, (tree, nodea, nodeb)-> iseven(nodea))))
```

and for the Petrov-Galerkin case

```@example filteringnodes2
using CompScienceMeshes # hide
using H2Trees # hide

mx = meshsphere(1.0, 0.1)
my = meshsphere(2.0, 0.1)

tree = TwoNTree(vertices(mx), vertices(my), 0.1)
testtree = H2Trees.testtree(tree)
trialtree = H2Trees.trialtree(tree)
println(collect(H2Trees.NodeFilterIterator(tree, 3, (testree, trialtree, testnode, trialnode)-> iseven(testnode))))
```

and

```@example filteringnodes3
using CompScienceMeshes # hide
using H2Trees # hide

mx = meshsphere(1.0, 0.1)
my = meshsphere(2.0, 0.1)

tree = TwoNTree(vertices(mx), vertices(my), 0.1)
testtree = H2Trees.testtree(tree)
trialtree = H2Trees.trialtree(tree)
println(collect(H2Trees.NodeFilterIterator(testtree, trialtree, 3, (testree, trialtree, testnode, trialnode)-> iseven(testnode))))
```
