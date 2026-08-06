"""
    DisaggregateTranslatePlan

Downward traversal plan for receiving translated moments.

The plan records which receiving nodes collect far-field translations and which
source nodes translate to them. It is paired with [`AggregatePlan`](@ref), which
produces the source moments.
"""
struct DisaggregateTranslatePlan{T} <: AbstractDisaggregationPlan
    translatingnodes::Vector{Dict{Int,Vector{Int}}} # Translating nodes
    receivingnodes_by_level::Vector{Vector{Int}}
    nodes::Vector{Vector{Int}} # Disaggregation nodes
    levels::UnitRange{Int} # Disaggregation levels
    isdisaggregationnode::Vector{Bool} # Does the node receive a moment direclty or via one of its parents ?
    rootoffset::Int # In case the tree is not rooted at 1
    tree::T
end

function DisaggregateTranslatePlan(
    translatingnodes::Vector{Dict{Int,Vector{Int}}},
    nodes::Vector{Vector{Int}},
    levels::UnitRange{Int},
    isdisaggregationnode::Vector{Bool},
    rootoffset::Int,
    tree,
)
    return _makedisaggregatetranslateplan(
        translatingnodes, nodes, levels, isdisaggregationnode, rootoffset, tree
    )
end

function plantranslationtrait(::DisaggregateTranslatePlan)
    return IsTranslatingPlan()
end

"""
    DisaggregateTranslatePlan(tree, TranslatingNodesIterator)
    DisaggregateTranslatePlan(testtree, trialtree, TranslatingNodesIterator)

Build a `DisaggregateTranslatePlan` from a translating-node iterator.

For a single tree, `TranslatingNodesIterator(node)` yields source nodes in
`tree` that translate to `node`.

For two trees, the plan is built on the first tree, while source nodes are
selected from the second. Block trees must be unwrapped explicitly with this
two-tree form.
"""
function DisaggregateTranslatePlan(tree, TranslatingNodesIterator)
    return DisaggregateTranslatePlan(tree, TranslatingNodesIterator, treetrait(tree))
end

function DisaggregateTranslatePlan(testtree, trialtree, TranslatingNodesIterator)
    return _builddisaggregatetranslateplan(testtree, trialtree, TranslatingNodesIterator)
end

function DisaggregateTranslatePlan(tree, TranslatingNodesIterator, ::isBlockTree)
    return throw(
        ArgumentError(
            "BlockTrees are not supported for DisaggregateTranslatePlan. Please specify which tree is used
            for the disaggregation.",
        ),
    )
end

function DisaggregateTranslatePlan(tree, TranslatingNodesIterator, ::AbstractTreeTrait)
    return _builddisaggregatetranslateplan(
        tree, _TranslatingFunctor(tree, TranslatingNodesIterator)
    )
end

function _builddisaggregatetranslateplan(testtree, trialtree, TranslatingNodesIterator)
    return _builddisaggregatetranslateplan(
        testtree,
        _TranslatingBlockTreeFunctor(testtree, trialtree, TranslatingNodesIterator),
    )
end

function mintranslationlevel(plan::DisaggregateTranslatePlan)
    return first(plan.levels)
end

function _makedisaggregatetranslateplan(
    translatingnodes, nodes, levels, isdisaggregationnode, rootoffset, tree
)
    return DisaggregateTranslatePlan(
        translatingnodes,
        _receivingnodes_by_level(translatingnodes),
        nodes,
        _validateddisaggregationlevels(levels),
        isdisaggregationnode,
        rootoffset,
        tree,
    )
end

function refreshreceivingnodes!(plan::DisaggregateTranslatePlan)
    resize!(plan.receivingnodes_by_level, length(plan.translatingnodes))
    for levelid in eachindex(plan.translatingnodes)
        plan.receivingnodes_by_level[levelid] = collect(
            keys(plan.translatingnodes[levelid])
        )
    end
    return plan
end

function _builddisaggregatetranslateplan(tree, TranslatingNodesIterator)
    return _buildtranslateplan(DisaggregateTranslateBuild(), tree, TranslatingNodesIterator)
end

function translatingnodes(plan::DisaggregateTranslatePlan, receivingnode::Int, level::Int)
    return plan[receivingnode, level]
end

function translatingnodes(plan::DisaggregateTranslatePlan)
    return plan.translatingnodes
end

function translatingnodes(plan::DisaggregateTranslatePlan, receivingnode::Int)
    return plan[receivingnode]
end

function receivingnodes(plan::DisaggregateTranslatePlan, level::Int)
    leveltolevelid(plan, level) < 1 && return Int[]
    leveltolevelid(plan, level) > length(plan.translatingnodes) && return Int[]
    # Cached as a `Vector` rather than a raw `Dict` `KeySet` -- see the matching
    # `AggregateTranslatePlan` method for why (feeds a `@tasks`-based parallel loop, whose
    # chunking needs an indexable collection).
    return plan.receivingnodes_by_level[leveltolevelid(plan, level)]
end

function Base.getindex(plan::DisaggregateTranslatePlan, receivingnode::Int, level::Int)
    level < mintranslationlevel(plan) && return Int[]
    tfnodes = plan.translatingnodes[leveltolevelid(plan, level)]

    if haskey(tfnodes, receivingnode)
        return tfnodes[receivingnode]
    else
        return Int[]
    end
end

function Base.getindex(plan::DisaggregateTranslatePlan, receivingnode::Int)
    level_ = level(plan.tree, receivingnode)
    return plan[receivingnode, level_]
end
