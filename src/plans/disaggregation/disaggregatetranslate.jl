"""
        DisaggregateTranslatePlan

Disaggregation translation plan for a single tree.

This is a translating plan used during downward traversal. The matching aggregation plan is an `AggregatePlan`.

Fields:

    - `translatingnodes`: Per level, a dictionary mapping receiving nodes (in `tree`) to
        source nodes, which do not have to be in `tree`.
    - `nodes`: Per level, nodes that must be visited during disaggregation.
    - `levels`: Contiguous disaggregation levels ordered from coarse to fine (root to leaves).
    - `isdisaggregationnode`: Boolean flag per node indicating whether the node is
        reached by the disaggregation traversal (either directly receiving translated data or
        via its ancestors).
    - `rootoffset`: Offset used to convert node ids to 1-based storage indices
    - `tree`: Tree for which the plan is built.
"""
struct DisaggregateTranslatePlan{T} <: AbstractDisaggregationPlan
    translatingnodes::Vector{Dict{Int,Vector{Int}}} # Translating nodes
    nodes::Vector{Vector{Int}} # Disaggregation nodes
    levels::UnitRange{Int} # Disaggregation levels
    isdisaggregationnode::Vector{Bool} # Does the node receive a moment direclty or via one of its parents ?
    rootoffset::Int # In case the tree is not rooted at 1
    tree::T
end

function plantranslationtrait(::DisaggregateTranslatePlan)
    return IsTranslatingPlan()
end

"""
        DisaggregateTranslatePlan(tree, TranslatingNodesIterator)
        DisaggregateTranslatePlan(testtree, trialtree, TranslatingNodesIterator)

Build a `DisaggregateTranslatePlan` from an iterator/functor that provides translating
source nodes for a receiving node.

Single-tree form:

  - `TranslatingNodesIterator(node)` yields source nodes in `tree` that translate to
    `node`.

Two-tree form:

  - `TranslatingNodesIterator(testnode)` is wrapped so the resulting plan is built on
    `testtree`, while source nodes are selected using `trialtree`.

For block trees, the tree used for disaggregation must be specified explicitly via the
two-tree constructor.
"""
function DisaggregateTranslatePlan(tree, TranslatingNodesIterator)
    return DisaggregateTranslatePlan(tree, TranslatingNodesIterator, treetrait(tree))
end

function DisaggregateTranslatePlan(testtree, trialtree, TranslatingNodesIterator)
    return _DisaggregateTranslatePlan(
        testtree,
        trialtree,
        TranslatingNodesIterator,
        treetrait(testtree),
        treetrait(trialtree),
    )
end

function DisaggregateTranslatePlan(tree, TranslatingNodesIterator, ::isBlockTree)
    return error(
        "BlockTrees are not supported for DisaggregateTranslatePlan. Please specify which tree is used
        for the disaggregation.",
    )
end

function DisaggregateTranslatePlan(tree, TranslatingNodesIterator, ::AbstractTreeTrait)
    return _DisaggregateTranslatePlan(
        tree, _TranslatingFunctor(tree, TranslatingNodesIterator)
    )
end

function _DisaggregateTranslatePlan(
    testtree, trialtree, TranslatingNodesIterator, ::AbstractTreeTrait, ::AbstractTreeTrait
)
    return _DisaggregateTranslatePlan(
        testtree,
        _TranslatingBlockTreeFunctor(testtree, trialtree, TranslatingNodesIterator);
    )
end

function mintranslationlevel(plan::DisaggregateTranslatePlan)
    return first(plan.levels)
end

function _DisaggregateTranslatePlan(tree, TranslatingNodesIterator)
    disaggregationlevels = zeros(Int, numberoflevels(tree))
    disaggregationnodes = Vector{Vector{Int}}(undef, length(disaggregationlevels))
    translatingnodes = Vector{Dict{Int,Vector{Int}}}(undef, length(disaggregationlevels))
    isdisaggregationnode = zeros(Bool, numberofnodes(tree))

    rootoffset = root(tree) - 1

    levels = collect(reverse(H2Trees.levels(tree)))

    lk = Threads.SpinLock()
    for level in levels
        leveldisaggregationnodes = Int[]
        leveltranslatingnodes = Dict{Int,Vector{Int}}()

        @threads for node in LevelIterator(tree, level)
            nodeindex = node - rootoffset
            # we only want to visit a node if it or one of its parents is receiving a far
            # field
            nodehastobevisited = false

            tfnodes = collect(Int, TranslatingNodesIterator(node))

            !isempty(tfnodes) && (nodehastobevisited = true)

            if !nodehastobevisited
                for parent in ParentUpwardsIterator(tree, node)
                    for _ in TranslatingNodesIterator(parent)
                        nodehastobevisited = true
                        break
                    end
                end
            end

            if nodehastobevisited
                lock(lk) do
                    push!(leveldisaggregationnodes, node)
                    if !isempty(tfnodes)
                        leveltranslatingnodes[nodeindex] = tfnodes
                    end
                end
                isdisaggregationnode[nodeindex] = true
            end
        end

        (isempty(leveldisaggregationnodes) && isempty(leveltranslatingnodes)) && continue
        disaggregationlevels[leveltolevelid(tree, level)] = level
        disaggregationnodes[leveltolevelid(tree, level)] = leveldisaggregationnodes
        translatingnodes[leveltolevelid(tree, level)] = leveltranslatingnodes
    end

    indicestodelete = Int[]
    for i in eachindex(disaggregationnodes)
        if !isassigned(disaggregationnodes, i)
            push!(indicestodelete, i)
        end
    end

    deleteat!(disaggregationlevels, indicestodelete)
    deleteat!(disaggregationnodes, indicestodelete)
    deleteat!(translatingnodes, indicestodelete)

    @assert disaggregationlevels == disaggregationlevels[begin]:disaggregationlevels[end]

    return DisaggregateTranslatePlan(
        translatingnodes,
        disaggregationnodes,
        disaggregationlevels[begin]:disaggregationlevels[end],
        isdisaggregationnode,
        rootoffset,
        tree,
    )
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
    return keys(plan.translatingnodes[leveltolevelid(plan, level)])
end

function Base.getindex(plan::DisaggregateTranslatePlan, receivingnode::Int, level::Int)
    level < H2Trees.mintranslationlevel(plan) && return Int[]
    tfnodes = plan.translatingnodes[leveltolevelid(plan, level)]

    if haskey(tfnodes, receivingnode)
        return tfnodes[receivingnode]
    else
        return Int[]
    end
end

function Base.getindex(plan::DisaggregateTranslatePlan, receivingnode::Int)
    level = H2Trees.level(plan.tree, receivingnode)
    return plan[receivingnode, level]
end
