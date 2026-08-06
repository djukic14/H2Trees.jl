"""
    _buildgalerkinplans(tree, aggregatenode, translatingnodesiterator, aggregationmode)

Construct the four Galerkin plans for `tree`:
`trialaggregationplan`, `testdisaggregationplan`, `testaggregationplan`, and
`trialdisaggregationplan`.

The returned [`PlanSet`](@ref) also contains `relevantlevels`, defined from the
minimum level at which aggregation and disaggregation are both meaningful up to
the deepest level of `tree`. Use `NamedTuple(plans)` when raw named-tuple
interoperability is required.

Internal: reached only through [`buildplans`](@ref).

Assumptions:

  - Aggregation and disaggregation are built on the same `tree` (`tree` is not a `BlockTree`).
"""
function _buildgalerkinplans(tree, aggregatenode, translatingnodesiterator, aggregationmode)
    trialaggregationplan, testdisaggregationplan = _galerkinplans(
        tree, aggregatenode, translatingnodesiterator, aggregationmode
    )
    testaggregationplan, trialdisaggregationplan = adjointplans(
        trialaggregationplan, testdisaggregationplan
    )

    minrelevantlevel = max(
        1,
        min(
            minaggregationlevel(trialaggregationplan),
            mindisaggregationlevel(testdisaggregationplan),
        ),
    )

    relevantlevels = minrelevantlevel:H2Trees.levels(tree)[end]

    return PlanSet(;
        trialaggregationplan=trialaggregationplan,
        testdisaggregationplan=testdisaggregationplan,
        testaggregationplan=testaggregationplan,
        trialdisaggregationplan=trialdisaggregationplan,
        relevantlevels=relevantlevels,
        tree=tree,
        family=GalerkinPlanFamily(),
    )
end

"""
    _galerkinplans(tree, aggregatenode, translatingnodesiterator, ::AggregateMode)

Build the forward Galerkin plan pair for aggregate-then-translate mode.

The upward side is an [`AggregatePlan`](@ref); far-field translations are stored
on the downward [`DisaggregateTranslatePlan`](@ref).
"""
function _galerkinplans(tree, aggregatenode, translatingnodesiterator, ::AggregateMode)
    trialaggregationplan = AggregatePlan(tree, aggregatenode(tree))
    testdisaggregationplan = DisaggregateTranslatePlan(tree, translatingnodesiterator(tree))

    return trialaggregationplan, testdisaggregationplan
end

"""
    _galerkinplans(tree, aggregatenode, translatingnodesiterator, ::AggregateTranslateMode)

Build the forward Galerkin plan pair for aggregate-translate mode.

The upward side is an [`AggregateTranslatePlan`](@ref); the matching downward
side is a non-translating [`DisaggregatePlan`](@ref).
"""
function _galerkinplans(
    tree, aggregatenode, translatingnodesiterator, ::AggregateTranslateMode
)
    trialaggregationplan = AggregateTranslatePlan(tree, translatingnodesiterator(tree))
    testdisaggregationplan = DisaggregatePlan(tree, aggregatenode(tree))

    return trialaggregationplan, testdisaggregationplan
end
