"""
    DisaggregatePlan

Downward traversal plan for already-translated moments.

The plan stores disaggregation nodes grouped by level, ordered from root toward
the leaves. It is paired with [`AggregateTranslatePlan`](@ref), which produces
the translated moments this plan propagates.
"""
struct DisaggregatePlan{T} <: AbstractDisaggregationPlan
    nodes::Vector{Vector{Int}} # Disaggregation nodes
    levels::UnitRange{Int} # Disaggregation levels
    storenode::Vector{Bool} # Does the node receive a moment directly?
    rootoffset::Int
    tree::T
end

"""
    DisaggregatePlan(tree, disaggregatenode)

Build a `DisaggregatePlan` from `disaggregatenode(node)`.

`disaggregatenode` marks nodes that store moments directly; nodes below a marked
ancestor are included in the disaggregation traversal even when they do not
store directly.

This constructor is defined for non-`BlockTree` trees. For block trees, select
the corresponding test or trial tree first.
"""
function DisaggregatePlan(tree, disaggregatenode)
    return DisaggregatePlan(tree, disaggregatenode, treetrait(tree))
end

function DisaggregatePlan(tree, disaggregatenode, ::isBlockTree)
    return throw(
        ArgumentError(
            "BlockTrees are not supported for DisaggregatePlan. Please specify which tree is used
            for the disaggregation.",
        ),
    )
end

function DisaggregatePlan(tree, disaggregatenode, ::AbstractTreeTrait)
    return _buildstorednodeplan(DisaggregateStoredNodeBuild(), tree, disaggregatenode)
end
