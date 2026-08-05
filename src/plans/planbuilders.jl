"""
    AbstractPlanFamily

Marker supertype for the plan family stored in a [`PlanSet`](@ref).
"""
abstract type AbstractPlanFamily end

"""
    GalerkinPlanFamily <: AbstractPlanFamily

Plan-family marker for single-tree Galerkin plans.
"""
struct GalerkinPlanFamily <: AbstractPlanFamily end

"""
    PetrovPlanFamily <: AbstractPlanFamily

Plan-family marker for two-tree Petrov-Galerkin plans.
"""
struct PetrovPlanFamily <: AbstractPlanFamily end

"""
    PlanBuilder(; aggregationmode=AggregateMode(), isnear=H2Trees.isnear(),
        translatingnodesiterator=WellSeparatedIterator(; isnear=isnear),
        aggregatenode=istranslatingnode(; TranslatingNodesIterator=translatingnodesiterator),
        checkadmissibility=false)

Construction settings for [`buildplans`](@ref).

`translatingnodesiterator` and `aggregatenode` both default from `isnear`, so overriding `isnear`
alone (rather than reconstructing all three by hand) keeps them consistent. Passing an explicit
`translatingnodesiterator`/`aggregatenode` overrides the chained default, e.g. to decouple the
near/far classification used for translation from the one used elsewhere.
"""
struct PlanBuilder{A,I,TNI,AN,C}
    aggregationmode::A
    isnear::I
    translatingnodesiterator::TNI
    aggregatenode::AN
    checkadmissibility::C
end

function PlanBuilder(;
    aggregationmode=AggregateMode(),
    isnear=isnear(),
    translatingnodesiterator=WellSeparatedIterator(; isnear=isnear),
    aggregatenode=istranslatingnode(; TranslatingNodesIterator=translatingnodesiterator),
    checkadmissibility=false,
)
    return PlanBuilder(
        aggregationmode, isnear, translatingnodesiterator, aggregatenode, checkadmissibility
    )
end

"""
    PlanSet(; trialaggregationplan, testdisaggregationplan, testaggregationplan,
        trialdisaggregationplan, relevantlevels, tree, mintranslationlevel=nothing,
        family)

Concrete container returned by [`buildplans`](@ref).

`PlanSet` behaves like the named tuple used by the older API: it supports
iteration, `keys`, `values`, `pairs`, symbol indexing, `NamedTuple(plans)`, and
`merge`. Prefer the accessor functions when writing new code.
"""
struct PlanSet{TrA,TeD,TeA,TrD,L,T,F}
    trialaggregationplan::TrA
    testdisaggregationplan::TeD
    testaggregationplan::TeA
    trialdisaggregationplan::TrD
    relevantlevels::L
    tree::T
    mintranslationlevel::Union{Int,Nothing}
    family::F
end

function PlanSet(;
    trialaggregationplan,
    testdisaggregationplan,
    testaggregationplan,
    trialdisaggregationplan,
    relevantlevels,
    tree,
    mintranslationlevel=nothing,
    family,
)
    return PlanSet(
        trialaggregationplan,
        testdisaggregationplan,
        testaggregationplan,
        trialdisaggregationplan,
        relevantlevels,
        tree,
        mintranslationlevel,
        family,
    )
end

"""
    trialaggregationplan(plans)
    testdisaggregationplan(plans)
    testaggregationplan(plans)
    trialdisaggregationplan(plans)

Return the four concrete plans in a [`PlanSet`](@ref).
"""
trialaggregationplan(plans::PlanSet) = plans.trialaggregationplan
testdisaggregationplan(plans::PlanSet) = plans.testdisaggregationplan
testaggregationplan(plans::PlanSet) = plans.testaggregationplan
trialdisaggregationplan(plans::PlanSet) = plans.trialdisaggregationplan

"""
    relevantlevels(plans)

Return the levels where all plans in `plans` are relevant to the traversal.
"""
relevantlevels(plans::PlanSet) = plans.relevantlevels

"""
    tree(plans)

Return the tree or block tree that `plans` were built from.
"""
tree(plans::PlanSet) = plans.tree

"""
    mintranslationlevel(plans)

Return the first level at which translations are relevant.

Galerkin plans do not store this separately, so the first relevant level is used.
"""
function mintranslationlevel(plans::PlanSet)
    return something(plans.mintranslationlevel, minimum(plans.relevantlevels))
end

"""
    planfamily(plans)

Return the family marker stored in `plans`.
"""
planfamily(plans::PlanSet) = plans.family

"""
    testtree(plans)
    trialtree(plans)

Return the test or trial tree associated with `plans`.

For Galerkin plans, both accessors resolve to the same underlying tree.
"""
testtree(plans::PlanSet) = testtree(plans.tree)
trialtree(plans::PlanSet) = trialtree(plans.tree)

"""
    isgalerkin(plans_or_family)
    ispetrov(plans_or_family)

Return whether a [`PlanSet`](@ref) or plan-family marker represents Galerkin or
Petrov-Galerkin plans.
"""
isgalerkin(plans::PlanSet) = isgalerkin(planfamily(plans))
ispetrov(plans::PlanSet) = ispetrov(planfamily(plans))
isgalerkin(::GalerkinPlanFamily) = true
isgalerkin(::PetrovPlanFamily) = false
ispetrov(::GalerkinPlanFamily) = false
ispetrov(::PetrovPlanFamily) = true

"""
    ownedtree(plan)

Return the tree owned by one concrete plan.
"""
ownedtree(plan::AbstractPlan) = tree(plan)

function _ownedandothertree(contexttree, plan::AbstractPlan)
    return _ownedandothertree(contexttree, plan, treetrait(contexttree))
end

function _ownedandothertree(contexttree, plan::AbstractPlan, ::AbstractTreeTrait)
    return (ownedtree(plan), contexttree)
end

function _ownedandothertree(contexttree, plan::AbstractPlan, ::isBlockTree)
    test = testtree(contexttree)
    trial = trialtree(contexttree)
    owned = ownedtree(plan)
    owned === test && return (test, trial)
    owned === trial && return (trial, test)
    return error(
        "Cannot resolve which side of the BlockTree this plan belongs to: the plan's tree is " *
        "identical to neither `testtree(tree)` nor `trialtree(tree)`. Pass the same BlockTree " *
        "the plans were built from.",
    )
end

"""
    receivingtree(contexttree, plan)
    translatingtree(contexttree, plan)

Resolve the receiving and translating trees for a translating plan.

For single-tree plans both trees are the same. For block-tree plans, `contexttree`
must be the block tree used to build the plans so the owned side can be matched
by object identity.
"""
function receivingtree(contexttree, plan::DisaggregateTranslatePlan)
    owned, _ = _ownedandothertree(contexttree, plan)
    return owned
end

function translatingtree(contexttree, plan::DisaggregateTranslatePlan)
    _, other = _ownedandothertree(contexttree, plan)
    return other
end

function receivingtree(contexttree, plan::AggregateTranslatePlan)
    _, other = _ownedandothertree(contexttree, plan)
    return other
end

function translatingtree(contexttree, plan::AggregateTranslatePlan)
    owned, _ = _ownedandothertree(contexttree, plan)
    return owned
end

function receivingtree(plans::PlanSet, plan::DisaggregateTranslatePlan)
    return receivingtree(tree(plans), plan)
end

function translatingtree(plans::PlanSet, plan::DisaggregateTranslatePlan)
    return translatingtree(tree(plans), plan)
end

function receivingtree(plans::PlanSet, plan::AggregateTranslatePlan)
    return receivingtree(tree(plans), plan)
end

function translatingtree(plans::PlanSet, plan::AggregateTranslatePlan)
    return translatingtree(tree(plans), plan)
end

function Base.propertynames(plans::PlanSet, private::Bool=false)
    names = (
        :trialaggregationplan,
        :testdisaggregationplan,
        :testaggregationplan,
        :trialdisaggregationplan,
        :relevantlevels,
        :tree,
        :family,
    )
    if ispetrov(plans)
        names = (names..., :mintranslationlevel)
    end
    private && return (names..., fieldnames(typeof(plans))...)
    return names
end

function _astuple(plans::PlanSet)
    if ispetrov(plans)
        return (
            testaggregationplan=plans.testaggregationplan,
            trialaggregationplan=plans.trialaggregationplan,
            testdisaggregationplan=plans.testdisaggregationplan,
            trialdisaggregationplan=plans.trialdisaggregationplan,
            relevantlevels=plans.relevantlevels,
            mintranslationlevel=plans.mintranslationlevel,
        )
    end

    return (
        trialaggregationplan=plans.trialaggregationplan,
        testdisaggregationplan=plans.testdisaggregationplan,
        testaggregationplan=plans.testaggregationplan,
        trialdisaggregationplan=plans.trialdisaggregationplan,
        relevantlevels=plans.relevantlevels,
    )
end

Base.keys(plans::PlanSet) = Base.keys(_astuple(plans))
Base.values(plans::PlanSet) = Base.values(_astuple(plans))
Base.pairs(plans::PlanSet) = Base.pairs(_astuple(plans))
Base.iterate(plans::PlanSet, state...) = iterate(_astuple(plans), state...)
Base.length(plans::PlanSet) = length(_astuple(plans))
Base.getindex(plans::PlanSet, key::Symbol) = getindex(_astuple(plans), key)
Base.NamedTuple(plans::PlanSet) = _astuple(plans)
Base.merge(nt::NamedTuple, plans::PlanSet) = merge(nt, _astuple(plans))
Base.merge(plans::PlanSet, nt::NamedTuple) = merge(_astuple(plans), nt)

"""
    buildplans(tree; builder=PlanBuilder())

Build Galerkin plans for a single tree or Petrov plans for a block tree.

`builder` controls the aggregation mode, near/far predicate, translating-node
iterator, stored-node predicate, and optional admissibility check. When
`builder.checkadmissibility` is true, the generated plan set is validated with
[`checkadmissibility`](@ref) before it is returned.
"""
function buildplans(tree; builder::PlanBuilder=PlanBuilder())
    plans = if treetrait(tree) isa isBlockTree
        _buildpetrovplans(
            tree,
            builder.aggregatenode,
            builder.translatingnodesiterator,
            builder.aggregationmode,
        )
    else
        _buildgalerkinplans(
            tree,
            builder.aggregatenode,
            builder.translatingnodesiterator,
            builder.aggregationmode,
        )
    end

    if builder.checkadmissibility
        report = checkadmissibility(tree, plans; isnear=builder.isnear)
        report.ok || error("built plans failed admissibility check: $report")
    end

    return plans
end
