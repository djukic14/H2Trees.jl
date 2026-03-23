# Galerkin Plans

In the Galerkin case, the test tree and trial tree are the same.
For distinct test and trial trees, see [Petrov Plans](petrovplans.md).

This page shows how to build Galerkin plans with [`galerkinplans`](@ref).
The single input space `X` is used for both test and trial roles, so one tree is sufficient.

```@example GalerkinPlans
using CompScienceMeshes 
using BEAST 
using H2Trees

m = meshicosphere(25, 1.0)
X = raviartthomas(m)
tree = TwoNTree(X, 0.0; minvalues=100)
```

## 1. Define the translating nodes iterator

```@example GalerkinPlans
tfiterator = H2Trees.TranslatingNodesIterator(;
    isnear=H2Trees.isnear(; additionalbufferboxes=1)
)

aggregatenode = H2Trees.istranslatingnode(; TranslatingNodesIterator=tfiterator)
```

- `tfiterator` defines how near and far interactions are separated.
- `aggregatenode` marks which nodes should be treated as translating nodes in the plan construction.

## 2. Construct plans in `AggregateMode()`

```@example GalerkinPlans
plans = H2Trees.galerkinplans(tree, aggregatenode, tfiterator, H2Trees.AggregateMode())
nothing #hide
```

In `AggregateMode()` the forward pair is:

- `trialaggregationplan`: `AggregatePlan`
- `testdisaggregationplan`: `DisaggregateTranslatePlan`

The transposed plans are built automatically:

- `testaggregationplan`: `AggregateTranslatePlan`
- `trialdisaggregationplan`: `DisaggregatePlan`

The `relevantlevels` field contains the levels where the operations are active.

## 3. Construct plans in `AggregateTranslateMode()`

```@example GalerkinPlans
plans = H2Trees.galerkinplans(tree, aggregatenode, tfiterator, H2Trees.AggregateTranslateMode())
nothing #hide
```

In `AggregateTranslateMode()` the forward pair is swapped:

- `trialaggregationplan`: `AggregateTranslatePlan`
- `testdisaggregationplan`: `DisaggregatePlan`

And the automatically built transposed plans become:

- `testaggregationplan`: `AggregatePlan`
- `trialdisaggregationplan`: `DisaggregateTranslatePlan`

The `relevantlevels` field again contains the levels where the operations are active.
