# Aggregate Translate Plan

```@raw html
<p style="text-align: center;">
    <img src="../../assets/aggregationtranslationtree.svg" alt="Aggregate Translate Tree" style="width: 60%;" />
</p>
```

[`AggregateTranslatePlan`](@ref) describes upward aggregation together with translations.
It specifies not only which aggregation nodes are visited, but also which receiving nodes
collect translated moments from those aggregation nodes.

This is a translating aggregation plan.
Its counterpart for a matrix-vector product is [`DisaggregatePlan`](disaggregateplan.md).

## Building an AggregateTranslatePlan

Conceptually, `AggregateTranslatePlan` is built from a tree and a translating-nodes iterator.
There are two supported forms:

- `AggregateTranslatePlan(tree, TranslatingNodesIterator)`
- `AggregateTranslatePlan(trialtree, testtree, TranslatingNodesIterator)`

Single-tree form:

- `TranslatingNodesIterator(node)` provides receiving nodes in the same tree.

Two-tree form:

- the plan is built on `trialtree`, so aggregation/source nodes in the plan are `trialtree` nodes,
- receiving nodes are selected with respect to `testtree`
- each value `receivingnodes[level][trialnode]` is a vector of translating nodes from `testtree`.

Nodes are visited if they translate directly, or if an ancestor translates.
This keeps the upward traversal complete for all translating contributions.

## Fields

An `AggregateTranslatePlan` stores:

- `receivingnodes`: per-level dictionary mapping receiving nodes to translating nodes.
    (two-tree form: keys are receiving nodes in `testtree`; values are translating nodes in `trialtree`).
- `nodes`: aggregation nodes grouped by level; these are nodes in the aggregation tree
    (two-tree form: `trialtree`).
- `levels`: aggregation levels, ordered leaves-to-root.
- `istranslatingnode`: booleans over aggregation-tree node indices indicating whether a node contributes as translating node
    (two-tree form: `trialtree`).
- `rootoffset`: offset used to convert nodes to compact 1-based indices.
- `tree`: the tree associated with the plan (two-tree form: `trialtree`).

## Example

The example below constructs Galerkin plans in `AggregateTranslateMode()` and extracts
the trial-side `AggregateTranslatePlan`.

```@example AggregateTranslateGalerkinPlans
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
plans = H2Trees.galerkinplans(tree, aggregatenode, tfiterator, H2Trees.AggregateTranslateMode())
aggregatetranslateplan = plans.trialaggregationplan
aggregatetranslateplan
```
