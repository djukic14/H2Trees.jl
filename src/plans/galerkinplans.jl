"""
    galerkinplans(tree, aggregatenode, translatingnodesiterator, aggregationmode)

Construct the four Galerkin plans for `tree`:
`trialaggregationplan`, `testdisaggregationplan`, `testaggregationplan`, and
`trialdisaggregationplan`.

The returned named tuple also contains `relevantlevels`, defined from the
minimum level at which aggregation and disaggregation are both meaningful up to
the deepest level of `tree`.

Assumptions:

  - Aggregation and disaggregation are built on the same `tree` (`tree` is not a `BlockTree`).
"""
function galerkinplans(tree, aggregatenode, translatingnodesiterator, aggregationmode)
    trialaggregationplan, testdisaggregationplan = _galerkinplans(
        tree, aggregatenode, translatingnodesiterator, aggregationmode
    )
    testaggregationplan, trialdisaggregationplan = H2Trees.adjointplans(
        trialaggregationplan, testdisaggregationplan
    )

    minrelevantlevel = max(
        1,
        min(
            H2Trees.minaggregationlevel(trialaggregationplan),
            H2Trees.mindisaggregationlevel(testdisaggregationplan),
        ),
    )

    relevantlevels = minrelevantlevel:H2Trees.levels(tree)[end]

    return (
        trialaggregationplan=trialaggregationplan,
        testdisaggregationplan=testdisaggregationplan,
        testaggregationplan=testaggregationplan,
        trialdisaggregationplan=trialdisaggregationplan,
        relevantlevels=relevantlevels,
    )
end

function _galerkinplans(tree, aggregatenode, translatingnodesiterator, ::AggregateMode)
    trialaggregationplan = AggregatePlan(tree, aggregatenode(tree))
    testdisaggregationplan = DisaggregateTranslatePlan(tree, translatingnodesiterator(tree))

    return trialaggregationplan, testdisaggregationplan
end

function _galerkinplans(
    tree, aggregatenode, translatingnodesiterator, ::AggregateTranslateMode
)
    trialaggregationplan = AggregateTranslatePlan(tree, translatingnodesiterator(tree))
    testdisaggregationplan = DisaggregatePlan(tree, aggregatenode(tree))

    return trialaggregationplan, testdisaggregationplan
end
