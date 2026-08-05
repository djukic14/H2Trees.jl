"""
    PetrovAggregationFunctor(aggregatenode, blocktree, notnodetree, nodetree)

Adapt a block-tree aggregation predicate to the one-tree `AggregatePlan` API.

`nodetree` is the tree being aggregated; `notnodetree` is the opposite side of
the Petrov block tree and is passed through so the predicate can evaluate
test/trial geometry together.
"""
struct PetrovAggregationFunctor{F,T,TNN,TN}
    aggregatenode::F
    blocktree::T
    notnodetree::TNN
    nodetree::TN
end

function (p::PetrovAggregationFunctor)(node::Int)
    return p.aggregatenode(p.blocktree)(p.notnodetree, p.nodetree, node)
end

"""
    PetrovDisaggregationFunctor(translatingnodesiterator, blocktree, notnodetree)

Adapt a block-tree translating-node iterator to the one-tree disaggregation
plan API.

The concrete disaggregation tree is supplied at call time; `notnodetree` is the
opposite side used to find translating nodes.
"""
struct PetrovDisaggregationFunctor{F,T,TNN}
    translatingnodesiterator::F
    blocktree::T
    notnodetree::TNN
end

function (p::PetrovDisaggregationFunctor)(nodetree, node::Int)
    return p.translatingnodesiterator(p.blocktree)(p.notnodetree, nodetree, node)
end

"""
    _buildpetrovplans(tree, aggregatenode, translatingnodesiterator, aggregationmode)

Construct the four Petrov-Galerkin plans for a block tree `tree`:
`testaggregationplan`, `trialaggregationplan`, `testdisaggregationplan`, and
`trialdisaggregationplan`.

The returned [`PlanSet`](@ref) also contains:

  - `relevantlevels`: levels where all generated plans are simultaneously relevant.
  - `mintranslationlevel`: first level where translations appear in both test and
    trial disaggregation plans.

Use `NamedTuple(plans)` when raw named-tuple interoperability is required.

Internal: reached only through [`buildplans`](@ref).
"""
function _buildpetrovplans(tree, aggregatenode, translatingnodesiterator, aggregationmode)
    trialtree_ = trialtree(tree)
    testtree_ = testtree(tree)

    trialaggregationplan, testdisaggregationplan = _petrovplans(
        tree,
        testtree_,
        trialtree_,
        aggregatenode,
        translatingnodesiterator,
        aggregationmode,
    )

    testaggregationplan, trialdisaggregationplan = adjointplans(
        trialaggregationplan, testdisaggregationplan
    )

    mintranslationlevel = min(
        mindisaggregationlevel(trialdisaggregationplan),
        mindisaggregationlevel(testdisaggregationplan),
    )

    minrelevantlevel = max(
        1,
        min(
            mintranslationlevel,
            minaggregationlevel(trialaggregationplan),
            minaggregationlevel(testaggregationplan),
        ),
    )

    lowerleaflevel = max(H2Trees.levels(trialtree_)[end], H2Trees.levels(testtree_)[end])
    relevantlevels = minrelevantlevel:lowerleaflevel

    return PlanSet(;
        testaggregationplan=testaggregationplan,
        trialaggregationplan=trialaggregationplan,
        testdisaggregationplan=testdisaggregationplan,
        trialdisaggregationplan=trialdisaggregationplan,
        relevantlevels=relevantlevels,
        mintranslationlevel=mintranslationlevel,
        tree=tree,
        family=PetrovPlanFamily(),
    )
end

"""
    _petrovplans(blocktree, testtree, trialtree, aggregatenode,
        translatingnodesiterator, ::AggregateMode)

Build the forward Petrov plan pair for aggregate-then-translate mode.

The trial tree is aggregated with an [`AggregatePlan`](@ref), and the test tree
receives translations through a [`DisaggregateTranslatePlan`](@ref).
"""
function _petrovplans(
    blocktree, testtree, trialtree, aggregatenode, translatingnodesiterator, ::AggregateMode
)
    trialaggregationplan = AggregatePlan(
        trialtree, PetrovAggregationFunctor(aggregatenode, blocktree, testtree, trialtree)
    )

    testdisaggregationplan = DisaggregateTranslatePlan(
        testtree,
        PetrovDisaggregationFunctor(translatingnodesiterator, blocktree, trialtree),
    )

    return trialaggregationplan, testdisaggregationplan
end

"""
    _petrovplans(blocktree, testtree, trialtree, aggregatenode,
        translatingnodesiterator, ::AggregateTranslateMode)

Build the forward Petrov plan pair for aggregate-translate mode.

The trial tree aggregates and translates with an [`AggregateTranslatePlan`](@ref);
the test tree is then handled by a non-translating [`DisaggregatePlan`](@ref).
"""
function _petrovplans(
    blocktree,
    testtree,
    trialtree,
    aggregatenode,
    translatingnodesiterator,
    ::AggregateTranslateMode,
)
    trialaggregationplan = AggregateTranslatePlan(
        trialtree,
        PetrovDisaggregationFunctor(translatingnodesiterator, blocktree, testtree),
    )

    testdisaggregationplan = DisaggregatePlan(
        testtree, PetrovAggregationFunctor(aggregatenode, blocktree, trialtree, testtree)
    )

    return trialaggregationplan, testdisaggregationplan
end
