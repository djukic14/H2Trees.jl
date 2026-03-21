"""
        AggregateTranslatePlan

Aggregation translation plan for a tree.

This is a translating plan used during upward traversal to compute and store moments that contribute to translation.
The matching disaggregation plan is a `DisaggregateTranslatePlan`.

Fields:

  - `receivingnodes`: Per level, a dictionary mapping receiving nodes (not necessarily in `tree`) to
    source aggregation node ids (in `tree`) whose moments translate to the receiver.
  - `nodes`: Per level, nodes that must be visited during aggregation (from leaves to root).
  - `levels`: Aggregation levels ordered from leaves to root.
  - `istranslatingnode`: Boolean flag per node index indicating whether the node
    contributes as a translating/source node.
  - `rootoffset`: Offset used to convert node ids to 1-based storage indices.
  - `tree`: Tree for which the plan is built.
"""
struct AggregateTranslatePlan{T} <: AbstractAggregationPlan
    receivingnodes::Vector{Dict{Int,Vector{Int}}} # receivingnodes[leveltolevelid(level)][disaggregationnode] = [translatingaggregationnodes]
    nodes::Vector{Vector{Int}} # Aggregation nodes
    levels::StepRange{Int,Int} # Aggregation levels
    istranslatingnode::Vector{Bool} # Does the aggregationnode translate its moment?
    rootoffset::Int # In case the tree is not rooted at 1
    tree::T
end

function plantranslationtrait(::AggregateTranslatePlan)
    return IsTranslatingPlan()
end

"""
        AggregateTranslatePlan(tree, TranslatingNodesIterator)
        AggregateTranslatePlan(testtree, trialtree, TranslatingNodesIterator)

Build an `AggregateTranslatePlan` from an iterator/functor that provides translating
target nodes for a source aggregation node.

Single-tree form:

  - `TranslatingNodesIterator(node)` yields receiving nodes in `tree` that receive
    translated information from `node`.

Two-tree form:

  - `TranslatingNodesIterator(testnode)` is wrapped so the resulting plan is built on
    `testtree`, while receiving nodes are selected using `trialtree`.

For block trees, the tree used for aggregation must be specified explicitly via the
two-tree constructor.
"""
function AggregateTranslatePlan(tree, TranslatingNodesIterator)
    return AggregateTranslatePlan(tree, TranslatingNodesIterator, treetrait(tree))
end

function AggregateTranslatePlan(testtree, trialtree, TranslatingNodesIterator)
    return _AggregateTranslatePlan(
        testtree,
        trialtree,
        TranslatingNodesIterator,
        treetrait(testtree),
        treetrait(trialtree),
    )
end

function AggregateTranslatePlan(tree, TranslatingNodesIterator, ::isBlockTree)
    return error(
        "BlockTrees are not supported for AggregateTranslatePlan. Please specify which tree is used
        for the disaggregation.",
    )
end

function AggregateTranslatePlan(tree, TranslatingNodesIterator, ::AbstractTreeTrait)
    return _AggregateTranslatePlan(
        tree, _TranslatingFunctor(tree, TranslatingNodesIterator)
    )
end

function _AggregateTranslatePlan(
    testtree, trialtree, TranslatingNodesIterator, ::AbstractTreeTrait, ::AbstractTreeTrait
)
    return _AggregateTranslatePlan(
        testtree,
        _TranslatingBlockTreeFunctor(testtree, trialtree, TranslatingNodesIterator);
    )
end

function mintranslationlevel(plan::AggregateTranslatePlan)
    return minlevel(plan)
end

function _AggregateTranslatePlan(tree, TranslatingNodesIterator)
    aggregationlevels = zeros(Int, numberoflevels(tree))
    aggregationnodes = Vector{Vector{Int}}(undef, length(aggregationlevels))
    receivingnodes = Vector{Dict{Int,Vector{Int}}}(undef, length(aggregationlevels))

    rootoffset = H2Trees.root(tree) - 1
    levels = collect(H2Trees.levels(tree))

    lk = Threads.SpinLock()
    for level in levels
        levelaggregationnodes = Int[]
        levelreceivingnodes = Dict{Int,Vector{Int}}()
        levelid = numberoflevels(tree) - level + minimumlevel(tree)

        @threads for node in LevelIterator(tree, level)
            nodehastobevisited = false

            tfnodes = collect(Int, TranslatingNodesIterator(node))

            !isempty(tfnodes) && (nodehastobevisited = true)

            if !nodehastobevisited
                for parent in ParentUpwardsIterator(tree, node)
                    for node in TranslatingNodesIterator(parent)
                        nodehastobevisited = true
                        break
                    end
                end
            end

            if nodehastobevisited
                lock(lk) do
                    push!(levelaggregationnodes, node)
                    for tfnode in tfnodes
                        if !haskey(levelreceivingnodes, tfnode)
                            levelreceivingnodes[tfnode] = [node]
                        else
                            push!(levelreceivingnodes[tfnode], node)
                        end
                    end
                end
            end

            (isempty(levelaggregationnodes) && isempty(levelreceivingnodes)) && continue
            aggregationlevels[levelid] = level
            aggregationnodes[levelid] = levelaggregationnodes
            receivingnodes[levelid] = levelreceivingnodes
        end
    end

    indicestodelete = Int[]
    for i in eachindex(aggregationnodes)
        if !isassigned(aggregationnodes, i)
            push!(indicestodelete, i)
        end
    end

    deleteat!(aggregationlevels, indicestodelete)
    deleteat!(aggregationnodes, indicestodelete)
    deleteat!(receivingnodes, indicestodelete)

    aggregationlevels = if !isempty(aggregationlevels)
        @assert aggregationlevels == maximum(aggregationlevels):-1:minimum(aggregationlevels)
        maximum(aggregationlevels):-1:minimum(aggregationlevels)
    end

    (isempty(aggregationnodes) || isempty(aggregationlevels)) &&
        error("Empty AggregatePlan not supported.")

    return AggregateTranslatePlan(
        receivingnodes,
        aggregationnodes,
        aggregationlevels,
        _computeistranslatingnodes(receivingnodes, tree),
        rootoffset,
        tree,
    )
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
    return keys(plan.receivingnodes[leveltolevelid(plan, level)])
end

function translatingnodes(plan::AggregateTranslatePlan, receivingnode::Int, level::Int)
    return plan.receivingnodes[leveltolevelid(plan, level)][receivingnode]
end

function Base.getindex(plan::AggregateTranslatePlan, receivingnode::Int, level::Int)
    level < H2Trees.mintranslationlevel(plan) && return Int[]
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
