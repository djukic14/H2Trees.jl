# Builder Workflow

Every tree is built with `buildtree` (or `buildforest`), configured by a builder object that
carries all construction settings. The builder's type selects the tree type; changing tree type
means passing a different builder, not calling a different function.

```julia
buildtree(points; builder=TwoNTreeBuilder())
buildtree(testpoints, trialpoints; builder=BlockTreeBuilder())   # -> BlockTree
buildforest(points, graph, weights; builder=MetisForestBuilder())
```

Direct constructors like `TwoNTree(points; builder=...)` are convenience wrappers around the
same builder; prefer `buildtree`/`buildforest` when the tree type should follow from the builder
alone.

## Builders

| Builder | Builds | Needs |
| --- | --- | --- |
| [`TwoNTreeBuilder`](@ref) | [`TwoNTree`](@ref) | — |
| [`BlockTreeBuilder`](@ref) | [`BlockTree`](@ref) | — |
| [`BoundingBallTreeBuilder`](@ref) | [`BoundingBallTree`](@ref) | — |
| [`KMeansTreeBuilder`](@ref) | `KMeansTree` | `ParallelKMeans.jl` |
| [`MetisTreeBuilder`](@ref) | `MetisTree` | `Metis.jl` |
| [`MetisForestBuilder`](@ref) | `MetisForest` | `Metis.jl` |
| [`SimpleHybridTreeBuilder`](@ref) | [`SimpleHybridTree`](@ref) | — (wraps an existing tree) |

The high-level builders have sane defaults: `TwoNTreeBuilder()`, `KMeansTreeBuilder()`,
`MetisTreeBuilder()`, and `MetisForestBuilder()` all work with no arguments. The generic
`BoundingBallTreeBuilder` is lower-level: pass an explicit splitter and split count, or use
`KMeansTreeBuilder`/`MetisTreeBuilder` for the packaged ball-tree strategies. See
[Tree Families](tree_families.md) for which one to pick.

## Example: same points, different builder

```@example builders1
using CompScienceMeshes # hide
using H2Trees

m = meshsphere(1.0, 0.1)
points = vertices(m)

adaptive = buildtree(points; builder=TwoNTreeBuilder(; minvalues=60))
uniform = buildtree(points; builder=TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0))
println(adaptive)
println(uniform)
```

Swapping in a clustering builder (with `ParallelKMeans.jl` loaded) needs no other code change:

```julia
using ParallelKMeans
clustered = buildtree(points; builder=KMeansTreeBuilder(; numberofclusters=8))
```

An unsupported `(input, builder)` combination throws an `ArgumentError` instead of silently
picking a default.
