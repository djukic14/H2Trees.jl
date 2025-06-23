function adjointplans(aggregationplan, disaggregationplan)
    return adjointaggregation(aggregationplan, disaggregationplan),
    adjointdisaggregation(aggregationplan, disaggregationplan)
end

# this plan aggregates and translates the tree of the disaggregation plan
function adjointaggregation(::AggregatePlan, disaggregationplan::DisaggregateTranslatePlan)
    return AggregateTranslatePlan(
        _inverttranslatingdict(translatingnodes(disaggregationplan)),
        reverse(disaggregationnodes(disaggregationplan)),
        reverse(disaggregationlevels(disaggregationplan)),
        rootoffset(disaggregationplan),
        tree(disaggregationplan),
    )
end

# this plan disaggregates the tree of the aggregation plan without translating
function adjointdisaggregation(
    aggregationplan::AggregatePlan, disaggregationplan::DisaggregateTranslatePlan
)
    aggregationtree = tree(aggregationplan)
    storenode = zeros(Bool, numberofnodes(aggregationtree))
    tfnodes = translatingnodes(disaggregationplan)

    lk = Threads.SpinLock()
    @threads for tfnodes in tfnodes
        for (_, nodes) in tfnodes
            for node in nodes
                nodeindex = node - rootoffset(aggregationplan)
                @lock lk storenode[nodeindex] = true
            end
        end
    end
    alevels = aggregationlevels(aggregationplan)
    @assert alevels == alevels[begin]:-1:alevels[end]

    return DisaggregatePlan(
        reverse(aggregationnodes(aggregationplan)),
        alevels[end]:alevels[begin],
        storenode,
        rootoffset(aggregationplan),
        aggregationtree,
    )
end

# this plan aggregates the tree of the disaggregation plan without translating
function adjointaggregation(::AggregateTranslatePlan, disaggregationplan::DisaggregatePlan)
    return AggregatePlan(
        reverse(disaggregationnodes(disaggregationplan)),
        reverse(disaggregationlevels(disaggregationplan)),
        storenode(disaggregationplan),
        rootoffset(disaggregationplan),
        tree(disaggregationplan),
    )
end

# this plan disaggregates and translates the tree of the aggregation plan
function adjointdisaggregation(aggregationplan::AggregateTranslatePlan, ::DisaggregatePlan)
    alevels = aggregationlevels(aggregationplan)
    @assert alevels == alevels[begin]:-1:alevels[end]

    dnodes = reverse(aggregationnodes(aggregationplan))
    offset = rootoffset(aggregationplan)
    isdisaggregationnode = zeros(Bool, numberofnodes(tree(aggregationplan)))

    for nodes in dnodes
        @threads for node in nodes
            isdisaggregationnode[node - offset] = true
        end
    end

    return DisaggregateTranslatePlan(
        _inverttranslatingdict(receivingnodes(aggregationplan)),
        dnodes,
        alevels[end]:alevels[begin],
        isdisaggregationnode,
        offset,
        tree(aggregationplan),
    )
end

function _inverttranslatingdict(rnodes)
    tfnodes = Vector{Dict{Int,Vector{Int}}}(undef, length(rnodes))

    @threads for i in eachindex(rnodes)
        nodes = rnodes[i]
        j = length(rnodes) - i + 1
        tfnodes[j] = Dict{Int,Vector{Int}}()
        for (key, values) in nodes
            for value in values
                if haskey(tfnodes[j], value)
                    push!(tfnodes[j][value], key)
                else
                    tfnodes[j][value] = [key]
                end
            end
        end
    end

    return tfnodes
end
