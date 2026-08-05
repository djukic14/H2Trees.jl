# Printing

Trees, forests, and diagnostic reports have two display forms: a compact one-line `show` (used
inline, e.g. in a `Vector` or `@show`) and a detailed multi-line `text/plain` summary (used at the
REPL / notebook top level).

```@example printing1
using CompScienceMeshes # hide
using H2Trees

m = meshsphere(1.0, 0.1)
tree = TwoNTree(vertices(m); builder=TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0))

sprint(show, tree)  # compact: type, nodes, leaves, level range, root
```

```@example printing1
tree  # text/plain: one row per level, with node/value/child counts and a tree-specific metric
```

The per-level metric is `halfsize` for box trees ([`TwoNTree`](@ref)) and `radius` for ball trees
(`KMeansTree`, [`BoundingBallTree`](@ref)).

## Forests

A [`Forest`](@ref) follows the same pattern: a compact one-liner (tree/node/leaf counts, level
range), or the detailed form listing every component tree — see [Forest](forest.md).

```julia
forest              # detailed: one summary line per component tree
sprint(show, forest) # compact
```

## Diagnostic reports

[`H2Trees.AdmissibilityReport`](@ref) has its own `text/plain` `show`, listing `ok`, error/warning
counts, and every finding with its geometric gap where relevant — see
[Admissibility Diagnostics](admissibility.md) for reading one.
