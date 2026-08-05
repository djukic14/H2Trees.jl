# Galerkin Plans

In the Galerkin case, the test tree and trial tree are the same.
For distinct test and trial trees, see [Petrov Plans](petrovplans.md).

This page shows how to build Galerkin plans with [`buildplans`](@ref).
The single input space `X` is used for both test and trial roles, so one tree is sufficient.
The result is a [`PlanSet`](@ref); see [Plans Overview](plans.md) for accessors
and named-tuple interoperability.

```@example GalerkinPlans
using CompScienceMeshes 
using BEAST 
using H2Trees

m = meshicosphere(25, 1.0)
X = raviartthomas(m)
tree = buildtree(X; builder=TwoNTreeBuilder(; minhalfsize=0.0, minvalues=100))
```

## 1. Construct plans in `AggregateMode()`

```@example GalerkinPlans
plans = H2Trees.buildplans(
    tree;
    builder=H2Trees.PlanBuilder(;
        isnear=H2Trees.isnear(; additionalbufferboxes=1),
        aggregationmode=H2Trees.AggregateMode(),
    ),
)
nothing #hide
```

`PlanBuilder`'s `translatingnodesiterator` and `aggregatenode` both default from `isnear`, so
overriding `isnear` alone (as above) keeps them consistent without reconstructing either by
hand. Pass `translatingnodesiterator`/`aggregatenode` explicitly only when they need to differ
from what `isnear` alone would produce.

In `AggregateMode()` the forward pair is:

- `trialaggregationplan`: `AggregatePlan`
- `testdisaggregationplan`: `DisaggregateTranslatePlan`

The transposed plans are built automatically:

- `testaggregationplan`: `AggregateTranslatePlan`
- `trialdisaggregationplan`: `DisaggregatePlan`

The `relevantlevels` field contains the levels where the operations are active.
Equivalent accessors such as `H2Trees.relevantlevels(plans)` and
`H2Trees.trialaggregationplan(plans)` are also available.

## 2. Construct plans in `AggregateTranslateMode()`

```@example GalerkinPlans
plans = H2Trees.buildplans(
    tree;
    builder=H2Trees.PlanBuilder(;
        isnear=H2Trees.isnear(; additionalbufferboxes=1),
        aggregationmode=H2Trees.AggregateTranslateMode(),
    ),
)
nothing #hide
```

In `AggregateTranslateMode()` the forward pair is swapped:

- `trialaggregationplan`: `AggregateTranslatePlan`
- `testdisaggregationplan`: `DisaggregatePlan`

And the automatically built transposed plans become:

- `testaggregationplan`: `AggregatePlan`
- `trialdisaggregationplan`: `DisaggregateTranslatePlan`

The `relevantlevels` field again contains the levels where the operations are active.
