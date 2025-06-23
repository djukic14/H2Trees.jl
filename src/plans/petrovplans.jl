struct PetrovAggregationFunctor{F,T,TNN,TN}
    aggregatenode::F
    blocktree::T
    notnodetree::TNN
    nodetree::TN
end

function (p::PetrovAggregationFunctor)(node::Int)
    return p.aggregatenode(p.blocktree)(p.notnodetree, p.nodetree, node)
end

struct PetrovDisaggregationFunctor{F,T,TNN}
    translatingnodesiterator::F
    blocktree::T
    notnodetree::TNN
end

function (p::PetrovDisaggregationFunctor)(nodetree, node::Int)
    return p.translatingnodesiterator(p.blocktree)(p.notnodetree, nodetree, node)
end

function petrovplans(tree, aggregatenode, translatingnodesiterator, aggregationmode)
    trialtree = H2Trees.trialtree(tree)
    testtree = H2Trees.testtree(tree)

    trialaggregationplan, testdisaggregationplan = _petrovplans(
        tree, testtree, trialtree, aggregatenode, translatingnodesiterator, aggregationmode
    )

    testaggregationplan, trialdisaggregationplan = H2Trees.adjointplans(
        trialaggregationplan, testdisaggregationplan
    )

    mintranslationlevel = min(
        H2Trees.mindisaggregationlevel(trialdisaggregationplan),
        H2Trees.mindisaggregationlevel(testdisaggregationplan),
    )

    minrelevantlevel = max(
        1,
        min(
            mintranslationlevel,
            H2Trees.minaggregationlevel(trialaggregationplan),
            H2Trees.minaggregationlevel(testaggregationplan),
        ),
    )

    lowerleaflevel = max(H2Trees.levels(trialtree)[end], H2Trees.levels(testtree)[end])
    relevantlevels = minrelevantlevel:lowerleaflevel

    return (
        testaggregationplan=testaggregationplan,
        trialaggregationplan=trialaggregationplan,
        testdisaggregationplan=testdisaggregationplan,
        trialdisaggregationplan=trialdisaggregationplan,
        relevantlevels=relevantlevels,
        mintranslationlevel=mintranslationlevel,
    )
end

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
