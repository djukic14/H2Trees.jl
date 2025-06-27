struct AggregateTranslatePlan{T} <: AbstractAggregationPlan
    receivingnodes::Vector{Dict{Int,Vector{Int}}}
    nodes::Vector{Vector{Int}}
    levels::StepRange{Int,Int}
    rootoffset::Int
    tree::T
end
#TODO: consider storing keys of receivingnodes dict for multithreaded looping

function plantranslationtrait(::AggregateTranslatePlan)
    return IsTranslatingPlan()
end

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
    translatingnodes = Vector{Dict{Int,Vector{Int}}}(undef, length(aggregationlevels))

    rootoffset = H2Trees.root(tree) - 1
    levels = collect(H2Trees.levels(tree))

    lk = Threads.SpinLock()
    for level in levels
        levelaggregationnodes = Int[]
        leveltranslatingnodes = Dict{Int,Vector{Int}}()
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
                        if !haskey(leveltranslatingnodes, tfnode)
                            leveltranslatingnodes[tfnode] = [node]
                        else
                            push!(leveltranslatingnodes[tfnode], node)
                        end
                    end
                end
            end

            (isempty(levelaggregationnodes) && isempty(leveltranslatingnodes)) && continue
            aggregationlevels[levelid] = level
            aggregationnodes[levelid] = levelaggregationnodes
            translatingnodes[levelid] = leveltranslatingnodes
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
    deleteat!(translatingnodes, indicestodelete)

    aggregationlevels = if !isempty(aggregationlevels)
        @assert aggregationlevels == maximum(aggregationlevels):-1:minimum(aggregationlevels)
        maximum(aggregationlevels):-1:minimum(aggregationlevels)
    end

    (isempty(aggregationnodes) || isempty(aggregationlevels)) &&
        error("Empty AggregatePlan not supported.")

    return AggregateTranslatePlan(
        translatingnodes, aggregationnodes, aggregationlevels, rootoffset, tree
    )
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
