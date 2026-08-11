# BlockTree and Petrov Trees

A [`BlockTree`](@ref) pairs a **test** tree and a **trial** tree, built from separate point sets.
It is the tree shape for Petrov-Galerkin problems, where the test and trial spaces differ. When
test and trial are the same space, build a single [`TwoNTree`](@ref) instead (the Galerkin case);
see [Builder Workflow](builders.md).

```@example blocktree1
using CompScienceMeshes # hide
using H2Trees

mx = meshsphere(1.0, 0.1)
my = meshsphere(2.0, 0.1)
translate!(my, [6.0, 0.0, 0.0])

tree = buildtree(
    vertices(mx),
    vertices(my);
    builder=BlockTreeBuilder(;
        test=TwoNTreeBuilder(; minhalfsize=0.1), trial=TwoNTreeBuilder(; minhalfsize=0.1)
    ),
)
```

`H2Trees.testtree(tree)` / `H2Trees.trialtree(tree)` give the two sides back. Translations run
from the trial tree (aggregation side) to the test tree (disaggregation side).

## Why both sides must agree on level scale

`test` and `trial` must share `minhalfsize` and the builder `root` id:  `BlockTreeBuilder`
validates this and throws an `ArgumentError` otherwise. This does **not** mean both sides have
the same geometric root box: their centers and root halfsizes may differ because the two point
sets may occupy different domains. The requirement is that comparable levels use the same box
size scale for near/far and translation purposes.

The two point sets are rarely the same physical size, though. To reconcile this without changing
`minhalfsize`, the smaller side is started at a *deeper* `minlevel`: its root box is subdivided
down until its halfsize matches what the larger side reaches at level 1. Both trees then line up
level-for-level from that point on, even though one has more levels overall.

```@example blocktree1
println("test starts at level ",  H2Trees.level(H2Trees.testtree(tree), H2Trees.root(H2Trees.testtree(tree))))
println("trial starts at level ", H2Trees.level(H2Trees.trialtree(tree), H2Trees.root(H2Trees.trialtree(tree))))
```

If you pass an explicit `minlevel` yourself, it must match the value this level-scale resolution
would have picked: a mismatch is an `ArgumentError`, not a silent override.

```@eval
using H2Trees
include(joinpath(pkgdir(H2Trees), "docs", "plotutils.jl"))
displayedcode(joinpath(pkgdir(H2Trees), "docs", "plots", "blocktree.jl"))
```

```@raw html
<object data="../assets/plots/blocktree.html" type="text/html"  style="width:100%; height:50vh;"> </object>
```

Test (pink) and trial (blue) boxes at the same level, from independently sized inputs.

See [Constructing Petrov Plans](plans/petrovplans.md) for building translation plans on a
`BlockTree`, and [Admissibility Diagnostics](admissibility.md) for validating the result.
