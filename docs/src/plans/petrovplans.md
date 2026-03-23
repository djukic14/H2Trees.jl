# Petrov Plans

This page shows how to build Petrov plans with [`petrovplans`](@ref).
In the Petrov case, the test tree and trial tree are generally different.
For identical test and trial trees, see [Galerkin Plans](galerkinplans.md).

Because test and trial spaces are different (`X` and `Y`), the construction uses distinct test and trial trees.

```@example PetrovPlans
using CompScienceMeshes
using BEAST
using H2Trees

mx = meshicosphere(25, 1.0)
my = meshicosphere(30, 2.0)
translate!(my, [0.0, 0.0, 5.0])

X = raviartthomas(mx)
Y = raviartthomas(my)

tree = TwoNTree(X, Y, 0.0; testminvalues=100, trialminvalues=100)
```

## 1. Define the translating nodes iterator

```@example PetrovPlans
tfiterator = H2Trees.TranslatingNodesIterator(;
    isnear=H2Trees.isnear(; additionalbufferboxes=1)
)

aggregatenode = H2Trees.istranslatingnode(; TranslatingNodesIterator=tfiterator)
```

- `tfiterator` defines how near and far interactions are separated.
- `aggregatenode` marks which nodes should be treated as translating nodes in the plan construction.

## 2. Construct plans in `AggregateMode()`

```@example PetrovPlans
plans = H2Trees.petrovplans(tree, aggregatenode, tfiterator, H2Trees.AggregateMode())
nothing #hide
```

In `AggregateMode()` the forward pair is:

- `trialaggregationplan`: `AggregatePlan`
- `testdisaggregationplan`: `DisaggregateTranslatePlan`

The transposed companion plans are built automatically:

- `testaggregationplan`: `AggregateTranslatePlan`
- `trialdisaggregationplan`: `DisaggregatePlan`

The `relevantlevels` field contains the levels where the operations are active.
The `mintranslationlevel` field gives the first level where translations occur.

## 3. Construct plans in `AggregateTranslateMode()`

```@example PetrovPlans
plans = H2Trees.petrovplans(tree, aggregatenode, tfiterator, H2Trees.AggregateTranslateMode())
nothing #hide
```

In `AggregateTranslateMode()` the forward pair is swapped:

- `trialaggregationplan`: `AggregateTranslatePlan`
- `testdisaggregationplan`: `DisaggregatePlan`

And the automatically built transposed plans become:

- `testaggregationplan`: `AggregatePlan`
- `trialdisaggregationplan`: `DisaggregateTranslatePlan`

The `relevantlevels` field contains the levels where the operations are active.
The `mintranslationlevel` field gives the first level where translations occur.
