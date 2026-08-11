# Node ID Semantics

Node ids in a [`TwoNTree`](@ref) are not arbitrary handles: their numeric order carries the
tree's physical layout. This page documents exactly what that layout guarantees, and (just as
importantly) what it does not.

## The layout

Bulk construction assigns final ids **level-major, Hilbert-ordered within each level**:

```text
Node IDs →

level Lmin:      [ nodes in Hilbert order ]
level Lmin + 1:  [ nodes in Hilbert order ]
level Lmin + 2:  [ nodes in Hilbert order ]
...
leaf level:      [ nodes in Hilbert order ]
                  ^                      ^
               first leaf             last leaf
```

Concretely:

  - the root keeps exactly the id it was asked for (`TwoNTreeBuilder(; root=...)`, default `1`),
    and all ids form the contiguous range `root : root + numberofnodes(tree) - 1`;
  - every level is one contiguous block of ids, so `nodesatlevel(tree, level)` is always a
    contiguous range;
  - adjacent levels touch: `last(nodesatlevel(tree, l)) + 1 == first(nodesatlevel(tree, l + 1))`;
  - within a level, increasing node id follows the Hilbert order of that level's occupied cells;
  - `firstchild`/`nextsibling` (and therefore [`H2Trees.children`](@ref)) also follow Hilbert order.

Ids are determined by geometry alone. Building the same point set in a different input order
produces the same node ids with the same geometry attached.

## Why level-major, and why Hilbert

Hierarchical methods do most of their work *between geometrically nearby boxes on the same
level*. Level-major storage makes level-wise kernels natural, and Hilbert ordering within the
level gives consecutive ids good spatial locality. Together they mean a fixed-size chunk of ids
tends to cover a compact region of space, so near-interaction data gets reused rather than
re-fetched.

## The curve itself

Before looking at trees, here is the plain Hilbert curve H2Trees orders by, in each supported
dimension. Marker colour runs along the traversal (dark → light), so the direction of travel is
visible and not just the shape of the path.

In **1D** the curve degenerates to the natural interval order: left to right at every
refinement, one orientation state, nothing to choose

```@eval
using H2Trees
include(joinpath(pkgdir(H2Trees), "docs", "plotutils.jl"))
displayedcode(joinpath(pkgdir(H2Trees), "docs", "plots", "hilbert_curve_1d.jl"))
```

```@raw html
<object data="../assets/plots/hilbert_curve_1d.html" type="text/html" style="width:100%; height:35vh;"> </object>
```

In **2D**, refined four times. Each level's path visits every cell exactly once, and consecutive
cells share a face: the property that makes consecutive ids spatially local:

```@eval
using H2Trees
include(joinpath(pkgdir(H2Trees), "docs", "plotutils.jl"))
displayedcode(joinpath(pkgdir(H2Trees), "docs", "plots", "hilbert_curve_2d.jl"))
```

```@raw html
<object data="../assets/plots/hilbert_curve_2d.html" type="text/html" style="width:100%; height:60vh;"> </object>
```

And in **3D** (drag to rotate). This is H2Trees' long-standing 3D convention, preserved exactly: 3D Hilbert curves are not unique, so the specific one is frozen rather than regenerated:

```@eval
using H2Trees
include(joinpath(pkgdir(H2Trees), "docs", "plotutils.jl"))
displayedcode(joinpath(pkgdir(H2Trees), "docs", "plots", "hilbert_curve_3d.jl"))
```

```@raw html
<object data="../assets/plots/hilbert_curve_3d.html" type="text/html" style="width:100%; height:40vh;"> </object>
```

## The curve through a real tree

The figures above are the abstract curve. What actually matters is that a built
[`TwoNTree`](@ref)'s node ids follow it: walking `nodesatlevel(tree, level)` (which is just
increasing node id) traces the curve through the node centers, with no sorting and no stored
Hilbert keys.

Each figure below builds a complete tree (one point per finest cell, so every box on every level
is occupied) and draws the polyline through the node centers in id order. **Use the slider to
step through the tree levels**; hover a marker to see its node id.

**1D** — ids simply run left to right:

```@eval
using H2Trees
include(joinpath(pkgdir(H2Trees), "docs", "plotutils.jl"))
displayedcode(joinpath(pkgdir(H2Trees), "docs", "plots", "hilbert_tree_1d.jl"))
```

```@raw html
<object data="../assets/plots/hilbert_tree_1d.html" type="text/html" style="width:100%; height:45vh;"> </object>
```

**2D** — the box grid is drawn faintly underneath, so the curve is visibly threading actual tree
cells rather than abstract grid points:

```@eval
using H2Trees
include(joinpath(pkgdir(H2Trees), "docs", "plotutils.jl"))
displayedcode(joinpath(pkgdir(H2Trees), "docs", "plots", "hilbert_tree_2d.jl"))
```

```@raw html
<object data="../assets/plots/hilbert_tree_2d.html" type="text/html" style="width:100%; height:65vh;"> </object>
```

**3D** — the deepest level here is the leaf level, so its path is exactly the order a
`Iterators.partition` over the leaf id range would chunk:

```@eval
using H2Trees
include(joinpath(pkgdir(H2Trees), "docs", "plotutils.jl"))
displayedcode(joinpath(pkgdir(H2Trees), "docs", "plots", "hilbert_tree_3d.jl"))
```

```@raw html
<object data="../assets/plots/hilbert_tree_3d.html" type="text/html" style="width:100%; height:60vh;"> </object>
```

!!! note
    These examples are deliberately *complete* grids. On a sparse tree only occupied boxes are
    stored, so the drawn path would skip cells: still Hilbert-ordered, but no longer
    face-adjacent step to step. See [Sparse levels](@ref "Sparse levels: locality, not strict adjacency")
    below.

### Deep, non-uniform trees

Real geometry is rarely a complete grid. The figures below use a dense core plus a sparse halo,
with `minvalues` stopping subdivision where the points thin out: so the tree runs much deeper
(9–15 levels here), only a fraction of each level's cells exist, and **leaves appear on many
different levels** rather than all at the bottom.

Stepping the slider shows what that costs and what survives. The path visibly *jumps* between
distant cells, because only occupied boxes are stored and the stored order is a subsequence of
the full curve: the caveat spelled out below. What still holds at every level is the ordering
itself: increasing node id follows the Hilbert order of the boxes that do exist. Each level's
title also reports how many of its nodes are leaves, which is the cross-level leaf distribution
that makes [`H2Trees.leaves`](@ref) non-contiguous for these trees.

**1D** (15 levels):

```@eval
using H2Trees
include(joinpath(pkgdir(H2Trees), "docs", "plotutils.jl"))
displayedcode(joinpath(pkgdir(H2Trees), "docs", "plots", "hilbert_tree_adaptive_1d.jl"))
```

```@raw html
<object data="../assets/plots/hilbert_tree_adaptive_1d.html" type="text/html" style="width:100%; height:45vh;"> </object>
```

**2D** (11 levels) — the deeper levels cover only the dense core, which is why the drawn boxes
shrink into one corner of the root:

```@eval
using H2Trees
include(joinpath(pkgdir(H2Trees), "docs", "plotutils.jl"))
displayedcode(joinpath(pkgdir(H2Trees), "docs", "plots", "hilbert_tree_adaptive_2d.jl"))
```

```@raw html
<object data="../assets/plots/hilbert_tree_adaptive_2d.html" type="text/html" style="width:100%; height:65vh;"> </object>
```

**3D** (9 levels):

```@eval
using H2Trees
include(joinpath(pkgdir(H2Trees), "docs", "plotutils.jl"))
displayedcode(joinpath(pkgdir(H2Trees), "docs", "plots", "hilbert_tree_adaptive_3d.jl"))
```

```@raw html
<object data="../assets/plots/hilbert_tree_adaptive_3d.html" type="text/html" style="width:100%; height:60vh;"> </object>
```

These are exactly the trees where the level-major layout earns its keep: see
[Adaptive trees](@ref "Adaptive trees: no globally contiguous leaf block") for the measured
effect on near-field chunking.

## Balanced trees: leaves are one contiguous chunkable range

When [`H2Trees.checkbalancedtree`](@ref) is `true`, every leaf sits on the deepest level. Since that
whole level is contiguous, the leaves are one contiguous Hilbert-ordered id range and can be
chunked directly:

```julia
leafids = leaves(tree)
@assert leafids == collect(first(leafids):last(leafids))

for chunk in Iterators.partition(first(leafids):last(leafids), chunksize)
    # geometrically local, Hilbert-ordered work
end
```

No permutation array, no stored Hilbert keys, no sorting at runtime: the ids themselves are the
ordering.

## What is *not* guaranteed

### Sparse levels: locality, not strict adjacency

Hilbert continuity gives face adjacency for consecutive cells of a **complete** grid. H2Trees
only stores occupied boxes, so a sparse level's stored order is a *subsequence* of the full
Hilbert traversal:

```text
full Hilbert grid:   A B C D E F G H I ...
stored boxes:        A B C       H I ...
resulting order:     A B C H I
```

`C` and `H` are consecutive ids but need not touch. The guarantee is:

> Increasing node ids follow the Hilbert order of the **occupied** boxes, and therefore have good
> spatial locality.

It is **not**:

> Every pair of consecutive node ids is geometrically adjacent.

### Adaptive trees: no globally contiguous leaf block

If leaves exist on several levels, it is impossible to have both "every level is one contiguous
id block" and "all leaves form one contiguous id block". H2Trees chooses **level-major**. So for
an adaptive tree each level is still contiguous and Hilbert-ordered, but `leaves(tree)` is spread
across several level blocks and is *not* a contiguous range. Only check for a contiguous leaf
range after confirming [`H2Trees.checkbalancedtree`](@ref).

Chunking an adaptive tree's leaves therefore needs one extra step — sort them by id first:

```julia
leafids = sort(leaves(tree))   # NOT a contiguous range, but level-grouped and Hilbert-ordered

for chunk in Iterators.partition(leafids, chunksize)
    # leaves of one level (mostly), in Hilbert order
end
```

This is worth doing rather than partitioning `leaves(tree)` directly, which comes back in
depth-first order. Near-field coupling is predominantly *same-level*, so grouping by level puts
same-scale boxes with heavily overlapping near fields into one chunk, while depth-first order
interleaves scales. Measured with `benchmark/locality.jl` (distinct near boxes per leaf, lower is
better), sorting by id needs **1.14x–1.63x fewer** distinct near boxes per chunk than depth-first
order across several adaptive geometries.

For balanced trees the two orders coincide, so the sort is harmless there and the plain
`first(leafids):last(leafids)` range above stays the cheapest option.

## Ids are only stable within one tree

The layout is a property of a constructed tree. Rebuilding with different builder settings
(different `minhalfsize`, `minvalues`, protrusion policy, …) changes the partition and therefore
the ids. Do not persist node ids across builds, and do not compare ids between two separately
built trees: compare geometry instead.

## Implementation notes

Hilbert mathematics lives in the internal `H2Trees.HilbertOrdering` submodule, which supports
dimensions 1, 2, and 3 and throws an `ArgumentError` outside that range. The 3D convention is the
long-standing H2Trees one (root sector order `0, 1, 3, 2, 6, 7, 5, 4`); 2D uses `0, 1, 3, 2`; 1D
is the natural interval order.

Renumbering happens exactly once, during bulk construction, before the tree and its
[`TreeIndex`](@ref) exist: so no caller can ever observe the temporary construction ids. It costs
`O(numberofnodes(tree))` and leaves no Hilbert state on the finished tree: normal tree use carries
no Hilbert overhead at all.
