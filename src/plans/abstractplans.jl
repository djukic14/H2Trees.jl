abstract type AbstractPlan end

abstract type AbstractAggregationPlan <: AbstractPlan end

abstract type AbstractDisaggregationPlan <: AbstractPlan end

abstract type AbstractPlanTranslationTrait end

struct IsTranslatingPlan <: AbstractPlanTranslationTrait end

struct IsNotTranslatingPlan <: AbstractPlanTranslationTrait end

function tree(plan::AbstractPlan)
    return plan.tree
end

function rootoffset(plan::AbstractPlan)
    return plan.rootoffset
end

function nodes(plan::AbstractPlan)
    return plan.nodes
end

function nodes(plan::AbstractPlan, level::Int)
    return nodes(plan)[leveltolevelid(plan, level)]
end

function levels(plan::AbstractPlan)
    return plan.levels
end

# AggregationPlan ##########################################################################

function minlevel(plan::AbstractAggregationPlan)
    return last(plan.levels)
end

function minaggregationlevel(plan::AbstractAggregationPlan)
    return minlevel(plan)
end

function aggregationnodes(plan::AbstractAggregationPlan, level::Int)
    return nodes(plan, level)
end

function aggregationnodes(plan::AbstractAggregationPlan)
    return nodes(plan)
end

function aggregationlevels(plan::AbstractAggregationPlan)
    return levels(plan)
end

function leveltolevelid(plan::AbstractAggregationPlan, level::Int)
    return length(plan.levels) - (level - minaggregationlevel(plan))
end

# DisaggregationPlan #######################################################################

function isdisaggregationnode(plan::AbstractDisaggregationPlan, node::Int)
    return plan.isdisaggregationnode[node - rootoffset(plan)]
end

function disaggregationlevels(plan::AbstractDisaggregationPlan)
    return levels(plan)
end

function disaggregationnodes(plan::AbstractDisaggregationPlan, level::Int)
    return nodes(plan, level)
end

function disaggregationnodes(plan::AbstractDisaggregationPlan)
    return nodes(plan)
end

function mindisaggregationlevel(plan::AbstractDisaggregationPlan)
    return first(plan.levels)
end

function leveltolevelid(plan::AbstractDisaggregationPlan, level::Int)
    return level - mindisaggregationlevel(plan) + 1
end

# PlanTranslationTrait #####################################################################

function istranslatingplan(plan::AbstractPlan)
    return plantranslationtrait(plan) isa IsTranslatingPlan
end

function plantranslationtrait(::AbstractPlan)
    return IsNotTranslatingPlan()
end

function storenode(plan)
    return storenode(plan, plantranslationtrait(plan))
end

function storenode(plan, ::IsNotTranslatingPlan)
    return plan.storenode
end

function storenode(plan, node::Int)
    return storenode(plan, node, plantranslationtrait(plan))
end

function storenode(plan, node::Int, ::IsNotTranslatingPlan)
    return storenode(plan)[node - rootoffset(plan)]
end

struct StoreNodeFunctor{A}
    plan::A
end

function (f::StoreNodeFunctor)(node::Int)
    return storenode(f.plan, node)
end

struct StoreNoNodeFunctor{A}
    plan::A
end

function (f::StoreNoNodeFunctor)(node::Int)
    return false
end

struct TranslatingFunctor{T,TF}
    tree::T
    translatingnodesiterator::TF
end

function (d::TranslatingFunctor)(node::Int)
    return d.translatingnodesiterator(d.tree, node)
end

struct TranslatingBlockTreeFunctor{TE,TR,TF}
    testtree::TE
    trialtree::TR
    translatingnodesiterator::TF
end

function (d::TranslatingBlockTreeFunctor)(testnode::Int)
    return d.translatingnodesiterator(d.trialtree, d.testtree, testnode)
end

function translatingplan(aggregationplan, disaggregationplan)
    return translatingplan(
        aggregationplan,
        disaggregationplan,
        plantranslationtrait(aggregationplan),
        plantranslationtrait(disaggregationplan),
    )
end

function translatingplan(
    ::A, disaggregationplan::D, ::IsNotTranslatingPlan, ::IsTranslatingPlan
) where {A<:AbstractAggregationPlan,D<:AbstractDisaggregationPlan}
    return disaggregationplan
end

function translatingplan(
    aggregationplan::A, ::D, ::IsTranslatingPlan, ::IsNotTranslatingPlan
) where {A<:AbstractAggregationPlan,D<:AbstractDisaggregationPlan}
    return aggregationplan
end

function translatingplan(
    ::A, ::D, ::Any, ::Any
) where {A<:AbstractAggregationPlan,D<:AbstractDisaggregationPlan}
    return error("Exactly one of the plans must be translating.")
end
