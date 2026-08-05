"""
    AggregateTranslatePlan

Upward traversal plan for translating/source moments.

The plan records which aggregation nodes must produce moments for far-field
translation and groups them by receiving node and level. The matching downward
non-translating plan is [`DisaggregatePlan`](@ref).
"""
struct AggregateTranslatePlan{T} <: AbstractAggregationPlan
    receivingnodes::Vector{Dict{Int,Vector{Int}}} # receivingnodes[leveltolevelid(level)][disaggregationnode] = [translatingaggregationnodes]
    receivingnodes_by_level::Vector{Vector{Int}}
    nodes::Vector{Vector{Int}} # Aggregation nodes
    levels::StepRange{Int,Int} # Aggregation levels
    istranslatingnode::Vector{Bool} # Does the aggregationnode translate its moment?
    rootoffset::Int # In case the tree is not rooted at 1
    tree::T
end

function AggregateTranslatePlan(
    receivingnodes::Vector{Dict{Int,Vector{Int}}},
    nodes::Vector{Vector{Int}},
    levels::StepRange{Int,Int},
    istranslatingnode::Vector{Bool},
    rootoffset::Int,
    tree,
)
    return _makeaggregatetranslateplan(
        receivingnodes, nodes, levels, istranslatingnode, rootoffset, tree
    )
end

function plantranslationtrait(::AggregateTranslatePlan)
    return IsTranslatingPlan()
end

"""
    AggregateTranslatePlan(tree, TranslatingNodesIterator)
    AggregateTranslatePlan(testtree, trialtree, TranslatingNodesIterator)

Build an `AggregateTranslatePlan` from a translating-node iterator.

For a single tree, `TranslatingNodesIterator(node)` yields the receiving nodes
that receive translated information from `node`.

For two trees, the plan is built on the first tree, while receiving nodes are
selected from the second. Block trees must be unwrapped explicitly with this
two-tree form.
"""
function AggregateTranslatePlan(tree, TranslatingNodesIterator)
    return AggregateTranslatePlan(tree, TranslatingNodesIterator, treetrait(tree))
end

function AggregateTranslatePlan(testtree, trialtree, TranslatingNodesIterator)
    return _buildaggregatetranslateplan(testtree, trialtree, TranslatingNodesIterator)
end

function AggregateTranslatePlan(tree, TranslatingNodesIterator, ::isBlockTree)
    return throw(
        ArgumentError(
            "BlockTrees are not supported for AggregateTranslatePlan. Please specify which tree is used
            for the disaggregation.",
        ),
    )
end

function AggregateTranslatePlan(tree, TranslatingNodesIterator, ::AbstractTreeTrait)
    return _buildaggregatetranslateplan(
        tree, _TranslatingFunctor(tree, TranslatingNodesIterator)
    )
end

function _buildaggregatetranslateplan(testtree, trialtree, TranslatingNodesIterator)
    return _buildaggregatetranslateplan(
        testtree,
        _TranslatingBlockTreeFunctor(testtree, trialtree, TranslatingNodesIterator),
    )
end

function mintranslationlevel(plan::AggregateTranslatePlan)
    return minlevel(plan)
end

function _receivingnodes_by_level(receivingnodes)
    return [collect(keys(levelreceivingnodes)) for levelreceivingnodes in receivingnodes]
end

function _makeaggregatetranslateplan(
    receivingnodes, nodes, levels, istranslatingnode, rootoffset, tree
)
    return AggregateTranslatePlan(
        receivingnodes,
        _receivingnodes_by_level(receivingnodes),
        nodes,
        _validatedaggregationlevels(levels),
        istranslatingnode,
        rootoffset,
        tree,
    )
end

function refreshreceivingnodes!(plan::AggregateTranslatePlan)
    resize!(plan.receivingnodes_by_level, length(plan.receivingnodes))
    for levelid in eachindex(plan.receivingnodes)
        plan.receivingnodes_by_level[levelid] = collect(keys(plan.receivingnodes[levelid]))
    end
    return plan
end

function _buildaggregatetranslateplan(tree, TranslatingNodesIterator)
    return _buildtranslateplan(AggregateTranslateBuild(), tree, TranslatingNodesIterator)
end

function _computeistranslatingnodes(receivingnodes, tree)
    istranslatingnodes = zeros(Bool, numberofnodes(tree))
    for receivingnodesdict in receivingnodes
        for (_, translatingnodes) in receivingnodesdict
            for translatingnode in translatingnodes
                istranslatingnodes[translatingnode - root(tree) + 1] = true
            end
        end
    end
    return istranslatingnodes
end

function receivingnodes(plan::AggregateTranslatePlan)
    return plan.receivingnodes
end

function receivingnodes(plan::AggregateTranslatePlan, level::Int)
    # Return the cached Vector, not Dict keys: threaded chunking downstream
    # requires an indexable collection.
    return plan.receivingnodes_by_level[leveltolevelid(plan, level)]
end

function translatingnodes(plan::AggregateTranslatePlan, receivingnode::Int, level::Int)
    return plan.receivingnodes[leveltolevelid(plan, level)][receivingnode]
end

function Base.getindex(plan::AggregateTranslatePlan, receivingnode::Int, level::Int)
    level < mintranslationlevel(plan) && return Int[]
    tfnodes = plan.receivingnodes[leveltolevelid(plan, level)]

    if haskey(tfnodes, receivingnode)
        return tfnodes[receivingnode]
    else
        return Int[]
    end
end

function istranslatingnode(plan::AggregateTranslatePlan, node::Int)
    return plan.istranslatingnode[node - plan.rootoffset]
end
