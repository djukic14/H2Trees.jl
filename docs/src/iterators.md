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

The [`WellSeparatedIterator`](@ref) is used to identify which translations should occur and which should not.
This determination is based on the concept of well-separated nodes.

!!! note
    Two nodes are considered well-separated if their parents are near each other and the nodes themselves are far apart.
    This assumes that child clusters are completely inside their parent clusters.

The [`WellSeparatedIterator`](@ref) can be configured using either an [`isnear`](@ref) or an [`iswellseparated`](@ref) function.

The following example demonstrates the usage of the [`WellSeparatedIterator`](@ref) in the Galerkin-case

```@example wellseparated1
using CompScienceMeshes # hide
using H2Trees # hide

m = meshsphere(1.0, 0.1)
tree = KMeansTree(vertices(m), 4; minvalues=60)

# Using the WellSeparatedIterator with the default isnear() function
println("First iterator\t", collect(H2Trees.WellSeparatedIterator(tree, 3)))

# Creating a functor without specifying a tree first
functor = H2Trees.WellSeparatedIterator()
iterator = functor(tree)
println("Second iterator\t", collect(iterator(tree, 3)))

# Specifying a custom iswellseparated criterion
functor = H2Trees.WellSeparatedIterator(; iswellseparated=(tree) -> (tree, nodea, nodeb) -> iseven(nodea))
iterator = functor(tree)
println("Third iterator\t", collect(iterator(tree, 3)))

# Specifying a custom isnear criterion
functor = H2Trees.WellSeparatedIterator(; isnear=(tree) -> (tree, nodea, nodeb) -> iseven(nodea))
iterator = functor(tree)
println("Fourth iterator\t", collect(iterator(tree, 3)))
```

When defining the [`isnear`](@ref) or [`iswellseparated`](@ref) criterion, it is necessary to provide a function that takes a tree as input and returns another function. This returned function is then used to evaluate the criterion.

At first glance, this may seem like an unnecessary layer of indirection. However, it actually provides a significant advantage: it enables precomputations that can be performed only once, when the criterion function is first created, rather than every time the iterator is called.

By allowing the initial function to perform any necessary precomputations and store the results, the returned function can then simply use these precomputed values to evaluate the criterion. This can significantly improve performance, especially when working with large trees or complex criteria.

For the Petrov-Galerkin case, we assume that translations occur from the `trialtree` to the `testtree`.
This scenario can be demonstrated with the following example

```@example wellseparated2
using CompScienceMeshes # hide
using H2Trees # hide

mx = meshsphere(1.0, 0.1)
my = meshsphere(2.0, 0.1)

tree = TwoNTree(vertices(mx), vertices(my), 0.1)

# Using the WellSeparatedIterator with the default isnear() function
println("First iterator\t", collect(H2Trees.WellSeparatedIterator(tree, 4)))

# Creating a functor without specifying a tree first
functor = H2Trees.WellSeparatedIterator()
iterator = functor(tree)
println("Second iterator\t", collect(iterator(H2Trees.testtree(tree), H2Trees.trialtree(tree), 4)))

# Specifying a custom iswellseparated criterion
functor = H2Trees.WellSeparatedIterator(;
    iswellseparated=(tree) -> (testtree, trialtree, testnode, trialnode) -> iseven(testnode)
)
iterator = functor(tree)
println("Third iterator\t", collect(iterator(H2Trees.testtree(tree), H2Trees.trialtree(tree), 4)))

# Specifying a custom isnear criterion
functor = H2Trees.WellSeparatedIterator(;
    isnear=(tree) -> (testtree, trialtree, testnode, trialnode) -> iseven(testnode)
)
iterator = functor(tree)
println("Fourth iterator\t", collect(iterator(H2Trees.testtree(tree), H2Trees.trialtree(tree), 4)))
```

In general, it is more efficient to use functors instead of functions in this context.

## Translating Nodes

The [`TranslatingNodesIterator`](@ref) is a wrapper for the [`WellSeparatedIterator`](@ref).

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
