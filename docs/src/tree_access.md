# Tree Access and Values

## Topology

| Function | Returns |
| --- | --- |
| `H2Trees.root(tree)` | the root node id |
| `H2Trees.parent(tree, node)` | the parent node id (`0` at the root) |
| `H2Trees.children(tree, node)` | iterator over direct children |
| `H2Trees.firstchild`/`H2Trees.nextsibling` | raw topology links `children` is built from |
| `H2Trees.isleaf(tree, node)` | whether `node` has no children |
| `H2Trees.leaves(tree, node=root(tree))` | leaf node ids below `node` |
| `H2Trees.nodesatlevel(tree, level)` | node ids at `level` |
| `H2Trees.depthfirstnodes(tree)` | cached depth-first node order |
| `H2Trees.level(tree, node)` / `H2Trees.levels(tree)` | a node's level / all levels present |

These, plus the [iterators](iterators.md), cover everything needed to walk a tree without
touching its internals directly.

### The cached index

Topology queries are backed by a cached [`TreeIndex`](@ref) ([`H2Trees.treeindex`](@ref)), built
once at construction. It is treated as read-only: if you mutate a tree's topology directly (rare
— normal use goes through the builder), call `H2Trees.rebuildtreeindex!(tree)` afterwards so
`leaves`, `nodesatlevel`, etc. stay correct.

## Values

A node's subtree stores a set of value ids (e.g. point or basis-function indices). Four functions
read them, trading off allocation for convenience:

| Function | Allocates | Use for |
| --- | --- | --- |
| `H2Trees.values(tree, node)` | yes, a fresh `Vector` | quick one-off inspection |
| [`H2Trees.appendvalues!`](@ref)`(out, tree, node)` | no — appends into caller-owned `out` | gathering across many nodes into one reused buffer |
| [`H2Trees.foreachvalue`](@ref)`(f, tree, node)` | no | running `f` over every value without materializing any vector |
| [`H2Trees.anyvalue`](@ref)`(f, tree, node)` | no, short-circuits | a yes/no predicate over a subtree |
| [`H2Trees.numberofvalues`](@ref)`(tree, node)` | no | just the count |

```@example treeaccess1
using CompScienceMeshes # hide
using H2Trees

m = meshsphere(1.0, 0.1)
tree = TwoNTree(vertices(m); builder=TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0))

# Count without allocating a vector:
H2Trees.numberofvalues(tree, H2Trees.root(tree))

# Gather into one reused buffer across several nodes:
out = Int[]
for node in H2Trees.LevelIterator(tree, 3)
    H2Trees.appendvalues!(out, tree, node)
end
length(out)

# Short-circuiting predicate:
H2Trees.anyvalue(v -> v == 1, tree, H2Trees.root(tree))
```

For anything beyond a quick check, prefer `appendvalues!`/`foreachvalue`/`anyvalue` over
`values(tree, node)` on internal nodes — `values` walks the same subtree but materializes a new
vector every call.
