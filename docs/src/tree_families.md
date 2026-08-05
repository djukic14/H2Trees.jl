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

```@example treefamilies1
using CompScienceMeshes, PlotlyJS # hide
using H2Trees
using ParallelKMeans

m = meshsphere(1.0, 0.1)
boxtree = TwoNTree(vertices(m); builder=TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0))
balltree = KMeansTree(vertices(m); builder=KMeansTreeBuilder(; numberofclusters=4, minvalues=60))

traces = [wireframe(skeleton(m, 1))]
for node in H2Trees.LevelIterator(boxtree, 4)
    push!(traces, H2Trees.tracecube(boxtree, node; mode="lines", line_color=:pink))
end
for node in H2Trees.LevelIterator(balltree, 3)
    push!(
        traces,
        H2Trees.traceball(
            balltree, node; colorscale=[[0, :blue], [1, :blue]], opacity=0.25, showscale=false
        ),
    )
end

p = PlotlyJS.plot(
    traces,
    Layout(;
        scene=attr(;
            xaxis=attr(; visible=false), yaxis=attr(; visible=false), zaxis=attr(; visible=false)
        ),
        showlegend=false,
    ),
)
savefig(p, "tree_families.html"); # hide
nothing #hide
```

```@raw html
<object data="../tree_families.html" type="text/html"  style="width:100%; height:50vh;"> </object>
```
