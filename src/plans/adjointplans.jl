"""
    adjointplans(aggregationplan, disaggregationplan)

Return the adjoint pair `(adjointaggregation, adjointdisaggregation)` associated
with the provided forward `aggregationplan` and `disaggregationplan`.

The forward pair must contain exactly one translating plan. The adjoint pair
keeps that invariant while swapping the translating role to the opposite
traversal direction.
"""
function adjointplans(aggregationplan, disaggregationplan)
    return adjointaggregation(aggregationplan, disaggregationplan),
    adjointdisaggregation(aggregationplan, disaggregationplan)
end

"""
    adjointaggregation(aggregationplan, disaggregationplan)

Build the aggregation-side plan used by the adjoint traversal.

For a forward `AggregatePlan`/`DisaggregateTranslatePlan` pair, the adjoint must
aggregate and translate on the forward disaggregation tree, so the translating
dictionary is inverted.
"""
function adjointaggregation(::AggregatePlan, disaggregationplan::DisaggregateTranslatePlan)
    receivingnodes = _inverttranslatingdict(translatingnodes(disaggregationplan))
    return _makeaggregatetranslateplan(
        receivingnodes,
        reverse(disaggregationnodes(disaggregationplan)),
        reverse(disaggregationlevels(disaggregationplan)),
        _computeistranslatingnodes(receivingnodes, tree(disaggregationplan)),
        rootoffset(disaggregationplan),
        tree(disaggregationplan),
    )
end

"""
    adjointdisaggregation(aggregationplan, disaggregationplan)

Build the disaggregation-side plan used by the adjoint traversal.

For a forward `AggregatePlan`/`DisaggregateTranslatePlan` pair, the adjoint
disaggregates the forward aggregation tree without translating.
"""
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
    _validatedaggregationlevels(alevels)

    return DisaggregatePlan(
        reverse(aggregationnodes(aggregationplan)),
        alevels[end]:alevels[begin],
        storenode,
        rootoffset(aggregationplan),
        aggregationtree,
    )
end

"""
    adjointaggregation(aggregationplan::AggregateTranslatePlan, disaggregationplan::DisaggregatePlan)

Build the non-translating aggregation plan for the adjoint of a forward
`AggregateTranslatePlan`/`DisaggregatePlan` pair.
"""
function adjointaggregation(::AggregateTranslatePlan, disaggregationplan::DisaggregatePlan)
    return AggregatePlan(
        reverse(disaggregationnodes(disaggregationplan)),
        reverse(disaggregationlevels(disaggregationplan)),
        storenode(disaggregationplan),
        rootoffset(disaggregationplan),
        tree(disaggregationplan),
    )
end

"""
    adjointdisaggregation(aggregationplan::AggregateTranslatePlan, disaggregationplan::DisaggregatePlan)

Build the translating disaggregation plan for the adjoint of a forward
`AggregateTranslatePlan`/`DisaggregatePlan` pair.

The forward aggregate-translation dictionary is inverted so adjoint receiving
nodes point back to the forward translating nodes.
"""
function adjointdisaggregation(aggregationplan::AggregateTranslatePlan, ::DisaggregatePlan)
    alevels = aggregationlevels(aggregationplan)
    _validatedaggregationlevels(alevels)

    dnodes = reverse(aggregationnodes(aggregationplan))
    offset = rootoffset(aggregationplan)
    isdisaggregationnode = zeros(Bool, numberofnodes(tree(aggregationplan)))

    for nodes in dnodes
        @threads for node in nodes
            isdisaggregationnode[node - offset] = true
        end
    end

    translatingnodes = _inverttranslatingdict(receivingnodes(aggregationplan))
    return _makedisaggregatetranslateplan(
        translatingnodes,
        dnodes,
        alevels[end]:alevels[begin],
        isdisaggregationnode,
        offset,
        tree(aggregationplan),
    )
end

"""
    _inverttranslatingdict(rnodes)

Invert a per-level translating-node dictionary.

The input maps `receivingnode => translatingnodes`; the result maps each
forward translating node to the forward receiving nodes that depend on it. Level
order is reversed because adjoint traversal runs in the opposite direction.
"""
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
