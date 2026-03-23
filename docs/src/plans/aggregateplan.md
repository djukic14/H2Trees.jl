# Aggregate Plan

```@raw html
<p style="text-align: center;">
    <img src="../../assets/aggregationtree.svg" alt="Aggregation Tree" style="width: 60%;" />
</p>
```

[`AggregatePlan`](@ref) describes the upward aggregation traversal on a tree.
It specifies:

- which nodes are visited per level during aggregation,
- which node moments must be stored for later phases, and

This is a non-translating aggregation plan.
Its counterpart for a matrix-vector product is [`DisaggregateTranslatePlan`](disaggregatetranslateplan.md).

## Building an AggregatePlan

Conceptually, `AggregatePlan` is built from a `tree` and a function `aggregatenode(node)`.

- If `aggregatenode(node)` is `true`, that node is marked for storage.
- If not, the node is still included in traversal when one of its ancestors is a storing node.
- Aggregation levels are organized from leaves toward the root.

This ensures aggregation paths remain complete for all moments that must be available later.

## Fields

An `AggregatePlan` stores:

- `nodes`: aggregation nodes grouped by level.
- `levels`: aggregation levels, ordered leaves-to-root.
- `storenode`: booleans indicating which node moments are stored.
- `rootoffset`: offset used to convert node ids to compact 1-based indices.
- `tree`: the tree associated with the plan.

## Example

The example below constructs Galerkin plans in `AggregateMode()` and extracts
the trial-side `AggregatePlan`.

```@example AggregateGalerkinPlans
using CompScienceMeshes #hide
using BEAST #hide
using H2Trees #hide
m = meshicosphere(25, 1.0) #hide
X = raviartthomas(m) #hide
tree = TwoNTree(X, 0.0; minvalues=100) #hide
tfiterator = H2Trees.TranslatingNodesIterator(; #hide
    isnear=H2Trees.isnear(; additionalbufferboxes=1) #hide
) #hide
aggregatenode = H2Trees.istranslatingnode(; TranslatingNodesIterator=tfiterator) #hide
plans = H2Trees.galerkinplans(tree, aggregatenode, tfiterator, H2Trees.AggregateMode())
aggregateplan = plans.trialaggregationplan
aggregateplan
```
