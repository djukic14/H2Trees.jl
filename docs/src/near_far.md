# Near and Far Predicates

`isnear(tree, nodea, nodeb)` decides whether two nodes are treated as near-field (direct
interaction) or far-field (translated). `isfar` is its negation. Everything that walks near/far
structure — [`NearNodeIterator`](@ref), [`FarNodeIterator`](@ref),
[`WellSeparatedIterator`](@ref), [`H2Trees.PlanBuilder`](@ref),
[`H2Trees.checkadmissibility`](@ref) — is built on the same predicate.

## Two shapes

A predicate is either **unresolved** — a factory taking the tree and returning the actual
comparison function — or **resolved** — the comparison function itself:

```@example nearfar1
using CompScienceMeshes # hide
using H2Trees

m = meshsphere(1.0, 0.1)
tree = TwoNTree(vertices(m); builder=TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0))

unresolved = H2Trees.isnear()       # factory: call with a tree to resolve
resolved = unresolved(tree)         # now callable on nodes directly
resolved(tree, 4, 5)
```

A plain user closure works the same way: `tree -> (tree, a, b) -> iseven(a)` is unresolved,
`(tree, a, b) -> iseven(a)` is already resolved. Functions receiving `isnear` (checkadmissibility,
`PlanBuilder`) accept either.

## One tree vs. two trees

For a single tree (Galerkin), the resolved shape is `(tree, testnode, trialnode) -> Bool`. For a
[`BlockTree`](@ref) (Petrov), it is `(testtree, trialtree, testnode, trialnode) -> Bool` instead —
passing the wrong arity raises a clear error rather than silently comparing the wrong nodes.

## The default geometric predicate

For box trees, `isnear` uses [`H2Trees.isneargap`](@ref): near when the axis-aligned gap between
the two boxes is at most `1.0 * min(halfsize)` ([`H2Trees.DEFAULTNEARGAPBOXES`](@ref)), rather
than a fixed grid neighborhood. This measures the *actual* geometric gap using each side's own
halfsize, so it stays correct even when the two sides of a `BlockTree` are independently built and
share no common grid. Ball trees use [`H2Trees.isnearradius`](@ref) instead, based on
center-distance and radius.

Avoid predicates that only work for one tree orientation — the gap/radius forms above are safe for
independently built trees; a predicate assuming a shared grid is not.

## Customizing

```@example nearfar1
plans = H2Trees.buildplans(
    tree; builder=H2Trees.PlanBuilder(; isnear=H2Trees.isnear(; additionalbufferboxes=1))
)
nothing #hide
```

`additionalbufferboxes` widens the near field by that many extra half-sizes; `PlanBuilder`'s
translating-nodes iterator and aggregation-node rule both derive from `isnear`, so overriding it
here keeps the whole plan consistent.
