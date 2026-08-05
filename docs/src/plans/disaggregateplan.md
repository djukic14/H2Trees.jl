# Disaggregate Plan

```@raw html
<p style="text-align: center;">
    <img src="../../assets/disaggregationtree.svg" alt="Disaggregate Tree" style="width: 60%;" />
</p>
```

[`DisaggregatePlan`](@ref) describes non-translating disaggregation traversal on a tree.
It specifies which nodes are visited while propagating moments from higher levels
toward lower levels.

This is a non-translating disaggregation plan.
Its counterpart for a matrix-vector product is [`AggregateTranslatePlan`](aggregatetranslateplan.md).

## Building a DisaggregatePlan

Conceptually, `DisaggregatePlan` is built from a `tree` and a function
`disaggregatenode(node)`.

- If `disaggregatenode(node)` is `true`, that node receives and stores a moment directly.
- If not, the node is still included in traversal when one of its ancestors is a storing node.
- Disaggregation levels are organized from root toward leaves.

This ensures disaggregation paths remain complete for all moments that must be
distributed to lower levels.

## Fields

A `DisaggregatePlan` stores:

- `nodes`: disaggregation nodes grouped by level.
- `levels`: contiguous disaggregation levels, ordered root-to-leaves.
- `storenode`: booleans indicating whether a node receives/stores directly.
- `rootoffset`: offset used to convert node ids to compact 1-based indices.
- `tree`: the tree associated with the plan.

## Example

The example below constructs Galerkin plans in `AggregateTranslateMode()` and
extracts the test-side `DisaggregatePlan`.

```@example DisaggregateGalerkinPlans
using CompScienceMeshes #hide
using BEAST #hide
using H2Trees #hide
m = meshicosphere(25, 1.0) #hide
X = raviartthomas(m) #hide
tree = buildtree(X; builder=TwoNTreeBuilder(; minhalfsize=0.0, minvalues=100)) #hide
plans = H2Trees.buildplans(
    tree;
    builder=H2Trees.PlanBuilder(;
        isnear=H2Trees.isnear(; additionalbufferboxes=1),
        aggregationmode=H2Trees.AggregateTranslateMode(),
    ),
)
disaggregateplan = plans.testdisaggregationplan
disaggregateplan
```
