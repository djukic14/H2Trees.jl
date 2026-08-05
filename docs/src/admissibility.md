# Admissibility Diagnostics

[`H2Trees.checkadmissibility`](@ref) validates an *assembled* set of translation plans against
the tree's own near/far geometry. It exists because a plan can be internally self-consistent —
every structural invariant holds — and still be geometrically wrong: two boxes that are actually
touching can end up scheduled as a far translation. `checkadmissibility` catches that class of
bug directly, by re-checking every scheduled pair against `isnear`.

Run it after changing a near/far predicate, a plan builder, or `BlockTree` construction — not on
every matvec; it is a diagnostic, not a hot-path check.

```@example admissibility1
using CompScienceMeshes # hide
using H2Trees

m = meshsphere(1.0, 0.1)
tree = TwoNTree(vertices(m); builder=TwoNTreeBuilder(; minhalfsize=0.0, minvalues=60))
plans = H2Trees.buildplans(tree; builder=H2Trees.PlanBuilder())

report = H2Trees.checkadmissibility(tree, plans; throw=false)
report.ok
```

## Reading a report

[`H2Trees.AdmissibilityReport`](@ref) has `ok` (`true` iff there are no `:error` findings — a
`:warning` alone does not clear it) and `findings`, a list of
[`H2Trees.AdmissibilityFinding`](@ref)s. Each finding has a `severity` (`:error`/`:warning`), a
`kind` (e.g. `:nearpairtranslated`, `:marginalgap`, `:coveragegap`), and a geometric `gap` when
applicable. [`H2Trees.haserrors`](@ref)/[`H2Trees.haswarnings`](@ref) check just one severity.

```@example admissibility1
report
```

## When it fails

Passing an `isnear` that disagrees with the plan's own predicate turns a legitimate translation
pair into a reported error:

```@example admissibility1
badreport = H2Trees.checkadmissibility(
    tree, plans; isnear=(tree, a, b) -> true, throw=false
)
badreport.ok
first(badreport.findings)
```

With `throw=true` (the default), errors raise instead of returning a report — use `throw=false`
to inspect findings programmatically. `coverage=true` (default) additionally checks that near
values plus plan-scheduled far values partition every leaf's targets exactly once; it is stronger
than recomputing the far field from [`FarNodeIterator`](@ref) but more expensive, so pass
`coverage=false` for a quick pass on large trees.

See [Near and Far Predicates](near_far.md) for the predicate shapes `isnear` accepts.
