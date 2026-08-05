# Petrov Plans

This page shows how to build Petrov plans with [`buildplans`](@ref).
In the Petrov case, the test tree and trial tree are generally different.
For identical test and trial trees, see [Galerkin Plans](galerkinplans.md).
The result is a [`PlanSet`](@ref); see [Plans Overview](plans.md) for accessors
and named-tuple interoperability.

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

tree = buildtree(
    X,
    Y;
    builder=BlockTreeBuilder(;
        test=TwoNTreeBuilder(; minhalfsize=0.0, minvalues=100),
        trial=TwoNTreeBuilder(; minhalfsize=0.0, minvalues=100),
    ),
)
```

## 1. Construct plans in `AggregateMode()`

```@example PetrovPlans
plans = H2Trees.buildplans(
    tree;
    builder=H2Trees.PlanBuilder(;
        isnear=H2Trees.isnear(; additionalbufferboxes=1),
        aggregationmode=H2Trees.AggregateMode(),
    ),
)
nothing #hide
```

`buildplans` dispatches on `tree`'s `treetrait`, so the same call builds Petrov plans here for
a `BlockTree`, and Galerkin plans for a plain tree (see [Galerkin Plans](galerkinplans.md)).
`PlanBuilder`'s `translatingnodesiterator` and `aggregatenode` both default from `isnear`, so
overriding `isnear` alone (as above) keeps them consistent without reconstructing either by hand.

In `AggregateMode()` the forward pair is:

- `trialaggregationplan`: `AggregatePlan`
- `testdisaggregationplan`: `DisaggregateTranslatePlan`

The transposed companion plans are built automatically:

- `testaggregationplan`: `AggregateTranslatePlan`
- `trialdisaggregationplan`: `DisaggregatePlan`

The `relevantlevels` field contains the levels where the operations are active.
The `mintranslationlevel` field gives the first level where translations occur.
Equivalent accessors such as `H2Trees.relevantlevels(plans)` and
`H2Trees.mintranslationlevel(plans)` are also available.

## 2. Construct plans in `AggregateTranslateMode()`

```@example PetrovPlans
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

The `relevantlevels` field contains the levels where the operations are active.
The `mintranslationlevel` field gives the first level where translations occur.
