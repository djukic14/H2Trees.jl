# Disaggregate Translate Plan

```@raw html
<p style="text-align: center;">
    <img src="../../assets/disaggregationtranslationtree.svg" alt="Disaggregate Translate Tree" style="width: 60%;" />
</p>
```

[`DisaggregateTranslatePlan`](@ref) describes disaggregation traversal together with
translations.
It specifies which disaggregation nodes are visited and, for each receiving node,
which translating nodes contribute to it.

This is a translating disaggregation plan.
Its counterpart for a matrix-vector product is [`AggregatePlan`](aggregateplan.md).

## Building a DisaggregateTranslatePlan

Conceptually, `DisaggregateTranslatePlan` is built from a tree and a
translating-nodes iterator.
There are two supported forms:

- `DisaggregateTranslatePlan(tree, TranslatingNodesIterator)`
- `DisaggregateTranslatePlan(testtree, trialtree, TranslatingNodesIterator)`

Single-tree form:

- `TranslatingNodesIterator(node)` provides translating nodes for receiving
    node `node` in the same tree.

Two-tree form:

- the plan is built on `testtree`, so receiving/disaggregation nodes in the plan are
    `testtree` nodes,
- translating nodes are selected with respect to `trialtree`,
- each value `translatingnodes[level][receivingnode]` is a vector of translating nodes
    from `trialtree`.

Nodes are visited if they receive translated data directly, or if an ancestor does.
This keeps the downward traversal complete for all translated contributions.

## Fields

A `DisaggregateTranslatePlan` stores:

- `translatingnodes`: per-level dictionary mapping receiving nodes to translating nodes.
    In two-tree form: keys are receiving nodes in `testtree`; values are translating nodes in `trialtree`.
- `nodes`: disaggregation nodes grouped by level; these are nodes in the
    disaggregation tree (two-tree form: `testtree`).
- `levels`: contiguous disaggregation levels, ordered root-to-leaves.
- `isdisaggregationnode`: booleans over disaggregation-tree node indices indicating
    whether a node is reached by disaggregation traversal.
- `rootoffset`: offset used to convert nodes to compact 1-based indices.
- `tree`: the tree associated with the plan (two-tree form: `testtree`).

## Example

The example below constructs Galerkin plans in `AggregateMode()` and extracts
the test-side `DisaggregateTranslatePlan`.

```@example DisaggregateTranslateGalerkinPlans
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
        aggregationmode=H2Trees.AggregateMode(),
    ),
)
disaggregatetranslateplan = plans.testdisaggregationplan
disaggregatetranslateplan
```
