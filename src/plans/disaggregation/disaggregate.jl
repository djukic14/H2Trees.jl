"""
        DisaggregatePlan

Non-translating disaggregation plan. The corresponding aggregation plan is a `AggregateTranslatePlan``

Fields:

    - `nodes`: disaggregation nodes grouped per level, ordered from lower to higher
        level (root to leaves).
    - `levels`: contiguous level range corresponding to `nodes`, with increasing
        order from root to leaves.
    - `storenode`: boolean marker per tree node indicating whether that node
        receives and stores a moment directly or if one of its ancestors does so.
    - `rootoffset`: offset used to map node ids to 1-based storage indices.
    - `tree`: tree on which disaggregation is defined.
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

    Build a `DisaggregatePlan` on `tree` using `disaggregatenode`.

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
    return error(
        "BlockTrees are not supported for DisaggregatePlan. Please specify which tree is used
        for the disaggregation.",
    )
end

function DisaggregatePlan(tree, disaggregatenode, ::AbstractTreeTrait)
    storenodesarray = zeros(Bool, numberofnodes(tree))
    root = H2Trees.root(tree)
    rootoffset = root - 1

    disaggregationlevels = zeros(Int, numberoflevels(tree))
    disaggregationnodes = Vector{Vector{Int}}(undef, length(disaggregationlevels))

    lk = Threads.SpinLock()
    levels = collect(reverse(H2Trees.levels(tree)))

    for level in levels
        leveldisaggregationnodes = Int[]

        @threads for node in LevelIterator(tree, level)
            nodeindex = node - root + 1

            if disaggregatenode(node)
                storenodesarray[nodeindex] = true
                @lock lk push!(leveldisaggregationnodes, node)
            end

            storenodesarray[nodeindex] && continue

            for parent in ParentUpwardsIterator(tree, node)
                if disaggregatenode(parent)
                    @lock lk push!(leveldisaggregationnodes, node)
                    break
                end
            end
        end

        isempty(leveldisaggregationnodes) && continue

        disaggregationlevels[leveltolevelid(tree, level)] = level
        disaggregationnodes[leveltolevelid(tree, level)] = leveldisaggregationnodes
    end

    indicestodelete = Int[]
    for i in eachindex(disaggregationnodes)
        if !isassigned(disaggregationnodes, i)
            push!(indicestodelete, i)
        end
    end

    deleteat!(disaggregationlevels, indicestodelete)
    deleteat!(disaggregationnodes, indicestodelete)

    @assert disaggregationlevels == disaggregationlevels[begin]:disaggregationlevels[end]

    return DisaggregatePlan(
        disaggregationnodes,
        disaggregationlevels[begin]:disaggregationlevels[end],
        storenodesarray,
        rootoffset,
        tree,
    )
end
