"""
    AbstractPlan

Common supertype for traversal plans.

Plans store tree nodes grouped by level and define the order in which a matrix
factorization traversal visits those nodes.
"""
abstract type AbstractPlan end

"""
    AbstractAggregationPlan <: AbstractPlan

Supertype for upward traversal plans.

Aggregation plans visit nodes from fine to coarse levels. Concrete plans may
either only aggregate stored moments or aggregate and translate in one pass.
"""
abstract type AbstractAggregationPlan <: AbstractPlan end

"""
    AbstractDisaggregationPlan <: AbstractPlan

Supertype for downward traversal plans.

Disaggregation plans visit nodes from coarse to fine levels. Concrete plans may
either only disaggregate stored moments or disaggregate and translate in one
pass.
"""
abstract type AbstractDisaggregationPlan <: AbstractPlan end

abstract type AbstractPlanTranslationTrait end

"""
    IsTranslatingPlan <: AbstractPlanTranslationTrait

Trait for plans that also carry translating-node data.
"""
struct IsTranslatingPlan <: AbstractPlanTranslationTrait end

"""
    IsNotTranslatingPlan <: AbstractPlanTranslationTrait

Trait for plans that only store traversal nodes.
"""
struct IsNotTranslatingPlan <: AbstractPlanTranslationTrait end

"""
    tree(plan)

Return the tree associated with `plan`.
"""
function tree(plan::AbstractPlan)
    return plan.tree
end

"""
    rootoffset(plan)

Return the node-number offset used by `plan` for indexing per-node storage.
"""
function rootoffset(plan::AbstractPlan)
    return plan.rootoffset
end

"""
    nodes(plan)
    nodes(plan, level)

Return all planned nodes, or the nodes stored for one tree level.
"""
function nodes(plan::AbstractPlan)
    return plan.nodes
end

function nodes(plan::AbstractPlan, level::Int)
    return nodes(plan)[leveltolevelid(plan, level)]
end

"""
    levels(plan)

Return the levels covered by `plan`, in the traversal order used by that plan.
"""
function levels(plan::AbstractPlan)
    return plan.levels
end

# AggregationPlan ##########################################################################

"""
    minlevel(plan::AbstractAggregationPlan)

Return the finest level covered by an aggregation plan.
"""
function minlevel(plan::AbstractAggregationPlan)
    return last(plan.levels)
end

"""
    minaggregationlevel(plan)

Return the finest aggregation level covered by `plan`.
"""
function minaggregationlevel(plan::AbstractAggregationPlan)
    return minlevel(plan)
end

"""
    aggregationnodes(plan)
    aggregationnodes(plan, level)

Return the aggregation nodes grouped by level, or the nodes at one level.
"""
function aggregationnodes(plan::AbstractAggregationPlan, level::Int)
    return nodes(plan, level)
end

function aggregationnodes(plan::AbstractAggregationPlan)
    return nodes(plan)
end

"""
    aggregationlevels(plan)

Return aggregation levels in traversal order, from fine to coarse.
"""
function aggregationlevels(plan::AbstractAggregationPlan)
    return levels(plan)
end

"""
    leveltolevelid(plan::AbstractAggregationPlan, level)

Map a tree level to the storage index used by an aggregation plan.
"""
function leveltolevelid(plan::AbstractAggregationPlan, level::Int)
    return length(plan.levels) - (level - minaggregationlevel(plan))
end

function _validatedaggregationlevels(levels)
    isempty(levels) && return levels

    expectedlevels = maximum(levels):-1:minimum(levels)
    if levels != expectedlevels
        throw(
            ArgumentError(
                "aggregation plan levels must be contiguous and descending, got $(collect(levels))",
            ),
        )
    end

    return expectedlevels
end

# DisaggregationPlan #######################################################################

"""
    isdisaggregationnode(plan, node)

Return whether `node` is marked as a stored disaggregation node.
"""
function isdisaggregationnode(plan::AbstractDisaggregationPlan, node::Int)
    return plan.isdisaggregationnode[node - rootoffset(plan)]
end

"""
    disaggregationlevels(plan)

Return disaggregation levels in traversal order, from coarse to fine.
"""
function disaggregationlevels(plan::AbstractDisaggregationPlan)
    return levels(plan)
end

"""
    disaggregationnodes(plan)
    disaggregationnodes(plan, level)

Return the disaggregation nodes grouped by level, or the nodes at one level.
"""
function disaggregationnodes(plan::AbstractDisaggregationPlan, level::Int)
    return nodes(plan, level)
end

function disaggregationnodes(plan::AbstractDisaggregationPlan)
    return nodes(plan)
end

"""
    mindisaggregationlevel(plan)

Return the coarsest level covered by a disaggregation plan.
"""
function mindisaggregationlevel(plan::AbstractDisaggregationPlan)
    return first(plan.levels)
end

"""
    leveltolevelid(plan::AbstractDisaggregationPlan, level)

Map a tree level to the storage index used by a disaggregation plan.
"""
function leveltolevelid(plan::AbstractDisaggregationPlan, level::Int)
    return level - mindisaggregationlevel(plan) + 1
end

function _validateddisaggregationlevels(levels)
    isempty(levels) && return levels

    expectedlevels = levels[begin]:levels[end]
    if levels != expectedlevels
        throw(
            ArgumentError(
                "disaggregation plan levels must be contiguous and ascending, got $(collect(levels))",
            ),
        )
    end

    return expectedlevels
end

# PlanTranslationTrait #####################################################################

"""
    istranslatingplan(plan)

Return `true` when `plan` stores translating-node data in addition to traversal
nodes.
"""
function istranslatingplan(plan::AbstractPlan)
    return plantranslationtrait(plan) isa IsTranslatingPlan
end

"""
    plantranslationtrait(plan)

Return the translation trait for `plan`.

Concrete translating plans override this to return [`IsTranslatingPlan`](@ref).
"""
function plantranslationtrait(::AbstractPlan)
    return IsNotTranslatingPlan()
end

"""
    storenode(plan)
    storenode(plan, node)

Return the stored-node mask for non-translating plans, or query whether `node`
is marked as stored.
"""
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

"""
    StoreNodeFunctor(plan)

Callable wrapper around `storenode(plan, node)`.
"""
struct StoreNodeFunctor{A}
    plan::A
end

function (f::StoreNodeFunctor)(node::Int)
    return storenode(f.plan, node)
end

"""
    StoreNoNodeFunctor(plan)

Callable predicate that never stores a node.
"""
struct StoreNoNodeFunctor{A}
    plan::A
end

function (f::StoreNoNodeFunctor)(node::Int)
    return false
end

struct _TranslatingFunctor{T,TF}
    tree::T
    translatingnodesiterator::TF
end

function (d::_TranslatingFunctor)(node::Int)
    return d.translatingnodesiterator(d.tree, node)
end

struct _TranslatingBlockTreeFunctor{TE,TR,TF}
    testtree::TE
    trialtree::TR
    translatingnodesiterator::TF
end

function (d::_TranslatingBlockTreeFunctor)(testnode::Int)
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

"""
    translatingplan(aggregationplan, disaggregationplan)

Return the one translating plan from a valid aggregation/disaggregation pair.

Exactly one side of a plan pair must carry translations: either
`AggregateTranslatePlan` with `DisaggregatePlan`, or `AggregatePlan` with
`DisaggregateTranslatePlan`.
"""
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
    return throw(ArgumentError("exactly one of the plans must be translating"))
end

"""
    AggregateAllNodesFunctor()

Predicate that marks every node as an aggregation node.
"""
struct AggregateAllNodesFunctor end

function (f::AggregateAllNodesFunctor)(node::Int)
    return true
end

function (f::AggregateAllNodesFunctor)(testtree, trialtree, trialnode::Int)
    return true
end

function (f::AggregateAllNodesFunctor)(tree)
    return f
end

function (f::AggregateAllNodesFunctor)(testtree, trialtree)
    return f
end

"""
    AggregateOnlyRootFunctor(tree)

Predicate that marks only `root(tree)` as an aggregation node.
"""
struct AggregateOnlyRootFunctor
    root::Int
    function AggregateOnlyRootFunctor(tree)
        return new(root(tree))
    end
end

function (f::AggregateOnlyRootFunctor)(node::Int)
    return node == f.root
end

function (f::AggregateOnlyRootFunctor)(tree)
    return AggregateOnlyRootFunctor(tree)
end

function (f::AggregateOnlyRootFunctor)(testtree, trialtree)
    return AggregateOnlyRootFunctor(testtree)
end

# Util functions ###########################################################################

function aggregateallnodes()
    return AggregateAllNodesFunctor()
end

function aggregaterootonly()
    return AggregateRootOnlyFunctor()
end

# Functors for utils #######################################################################
struct AggregateRootOnlyFunctor end

function (a::AggregateRootOnlyFunctor)(tree)
    return a(tree, H2Trees.treetrait(tree))
end

function (a::AggregateRootOnlyFunctor)(tree, ::Any)
    return AggregateRootNotBlockTreeFunctor(tree)
end

function (a::AggregateRootOnlyFunctor)(tree, ::isBlockTree)
    return AggregateRootBlockTreeFunctor(tree)
end

const AggregateRootFunctor = AggregateRootOnlyFunctor

struct AggregateRootNotBlockTreeFunctor{T}
    tree::T
end

function (f::AggregateRootNotBlockTreeFunctor)(node::Int)
    return root(f.tree) == node
end

struct AggregateRootBlockTreeFunctor{T}
    tree::T
end

function (::AggregateRootBlockTreeFunctor)(testtree, trialtree, trialnode::Int)
    return root(trialtree) == trialnode
end
