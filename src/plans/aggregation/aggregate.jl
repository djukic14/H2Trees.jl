"""
        AggregatePlan

Aggregation plan for a tree.

This plan specifies which nodes are visited during upward traversal and which node
moments must be stored for later use by translation/disaggregation phases.

This is not a translating plan and the corresponding translating plan is a `DisaggregateTranslatePlan`.

Fields:

    - `nodes`: Aggregation nodes grouped per level sorted from leave nodes to root.
    - `levels`: Aggregation levels (leaves to root).
    - `storenode`: Boolean flag per node index indicating whether that node's moment is
        stored.
    - `rootoffset`: Offset used to convert node ids to compact 1-based storage indices
    - `tree`: Tree for which the plan is built.
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

A node is marked for storage when `aggregatenode(node)` is `true`. Nodes that do not
store moments are still included in traversal if any ancestor satisfies
`aggregatenode`, ensuring paths toward storing nodes are aggregated correctly.

Block trees are not supported by this constructor; specify the aggregation tree
explicitly in a non-block-tree context.
"""
function AggregatePlan(tree, aggregatenode)
    return AggregatePlan(tree, aggregatenode, treetrait(tree))
end

function AggregatePlan(tree, aggregatenode, ::isBlockTree)
    return error(
        "BlockTrees are not supported for AggregatePlan. Please specify which tree is used
        for the aggregation."
    )
end

function AggregatePlan(tree, aggregatenode, ::AbstractTreeTrait)
    storenodesarray = zeros(Bool, numberofnodes(tree))
    root = H2Trees.root(tree)

    aggregationlevels = zeros(Int, numberoflevels(tree))
    aggregationnodes = Vector{Vector{Int}}(undef, length(aggregationlevels))

    lk = Threads.SpinLock()
    levels = collect(H2Trees.levels(tree))

    for level in levels
        levelaggregationnodes = Int[]
        levelid = numberoflevels(tree) - level + minimumlevel(tree)

        # @threads
        for node in LevelIterator(tree, level)
            nodeindex = node - root + 1

            if aggregatenode(node)
                storenodesarray[nodeindex] = true
                @lock lk push!(levelaggregationnodes, node)
            end

            storenodesarray[nodeindex] && continue

            for parent in ParentUpwardsIterator(tree, node)
                if aggregatenode(parent)
                    @lock lk push!(levelaggregationnodes, node)
                    break
                end
            end
        end

        isempty(levelaggregationnodes) && continue

        aggregationlevels[levelid] = level
        aggregationnodes[levelid] = levelaggregationnodes
    end

    indicestodelete = Int[]
    for i in eachindex(aggregationnodes)
        if !isassigned(aggregationnodes, i)
            push!(indicestodelete, i)
        end
    end

    deleteat!(aggregationlevels, indicestodelete)
    deleteat!(aggregationnodes, indicestodelete)

    aggregationlevels = if !isempty(aggregationlevels)
        @assert aggregationlevels == maximum(aggregationlevels):-1:minimum(aggregationlevels)
        maximum(aggregationlevels):-1:minimum(aggregationlevels)
    end

    (isempty(aggregationnodes) || isempty(aggregationlevels)) &&
        error("Empty AggregatePlan not supported.")

    return AggregatePlan(
        aggregationnodes, aggregationlevels, storenodesarray, root - 1, tree
    )
end
