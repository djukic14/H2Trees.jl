# Tree Families

| Tree | Geometry | Builder | Needs | Deterministic | Use when |
| --- | --- | --- | --- | --- | --- |
| [`TwoNTree`](@ref) | axis-aligned boxes, $2^n$ children | [`TwoNTreeBuilder`](@ref) | — | yes | grid-like geometry; protrusion control matters (e.g. BEAST basis functions) |
| [`BoundingBallTree`](@ref) | bounding balls | [`BoundingBallTreeBuilder`](@ref) | — | depends on `splitter` | generic clustering with a custom splitter |
| `KMeansTree` | balls from k-means clusters | [`KMeansTreeBuilder`](@ref) | `ParallelKMeans.jl` | yes (seeded `rng` by default) | irregular point clouds without grid structure |
| `MetisTree` | balls from graph partitions | [`MetisTreeBuilder`](@ref) | `Metis.jl` | yes | mesh/graph connectivity should drive clustering |
| `MetisForest` | one `MetisTree` per component | [`MetisForestBuilder`](@ref) | `Metis.jl` | yes | disconnected geometry (e.g. multiple bodies) |
| [`SimpleHybridTree`](@ref) | wraps a `TwoNTree`, splits upper/lower levels | [`SimpleHybridTreeBuilder`](@ref) | — | yes | different treatment above/below a chosen level |
| [`BlockTree`](@ref) | a pair of trees (test + trial) | [`BlockTreeBuilder`](@ref) | — | yes | Petrov-Galerkin, test space ≠ trial space — see [BlockTree and Petrov Trees](blocktree_petrov.md) |

All builders are constructed through [`buildtree`](@ref)/[`buildforest`](@ref) — see
[Builder Workflow](builders.md).

## Box vs. ball clustering

The same points, split into boxes ([`TwoNTree`](@ref), pink) and balls (`KMeansTree`, blue) at
the same level:

```@eval
using H2Trees
include(joinpath(pkgdir(H2Trees), "docs", "plotutils.jl"))
displayedcode(joinpath(pkgdir(H2Trees), "docs", "plots", "tree_families.jl"))
```

```@raw html
<object data="../assets/plots/tree_families.html" type="text/html"  style="width:100%; height:50vh;"> </object>
```
