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
