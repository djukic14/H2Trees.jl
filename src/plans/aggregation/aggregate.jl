"""
    AggregatePlan

Upward traversal plan for one tree.

The plan stores aggregation nodes grouped by level, ordered from leaves toward
the root. `storenode` marks the nodes whose moments must persist after
aggregation for later translation. The matching downward translation plan is
[`DisaggregateTranslatePlan`](@ref).
"""
struct AggregatePlan{T} <: AbstractAggregationPlan
    nodes::Vector{Vector{Int}} # aggregation nodes sorted by level
    levels::StepRange{Int,Int} # aggregation levels
    storenode::Vector{Bool} # Does the moment of the box need to be stored?
    rootoffset::Int # In case the tree is not rooted at 1
    tree::T
end

"""
    AggregatePlan(tree, aggregatenode)

Build an `AggregatePlan` from `aggregatenode(node)`.

A node is stored when `aggregatenode(node)` is `true`. Non-stored nodes are
still visited when they lie below a stored ancestor, so all required upward paths
are included.

Block trees are not supported by this constructor; specify the aggregation tree
explicitly in a non-block-tree context.
"""
function AggregatePlan(tree, aggregatenode)
    return AggregatePlan(tree, aggregatenode, treetrait(tree))
end

function AggregatePlan(tree, aggregatenode, ::isBlockTree)
    return throw(
        ArgumentError(
            "BlockTrees are not supported for AggregatePlan. Please specify which tree is used
            for the aggregation.",
        ),
    )
end

function AggregatePlan(tree, aggregatenode, ::AbstractTreeTrait)
    return _buildstorednodeplan(AggregateStoredNodeBuild(), tree, aggregatenode)
end
