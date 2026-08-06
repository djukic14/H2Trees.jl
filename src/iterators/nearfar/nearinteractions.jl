
"""
    nearinteractions(tree; kwargs...)

Compute near-interaction index lists for `tree`.

For a single tree, returns `(values, nearvalues)`: each `values[i]` is paired
with `nearvalues[i]`. Self-interactions are included by default.

Set `extractselfvalues=true` to return
`(selfvalues, values, nearvalues)`, where self-interactions are separated from
non-self near interactions. Pass `isnear` to choose a custom near predicate.
"""
function nearinteractions(tree; kwargs...)
    return nearinteractions(tree, treetrait(tree); kwargs...)
end

"""
    nearinteractions(tree, ::isBlockTree; kwargs...)

Compute near interactions for the two sides of a block tree.

This delegates to `nearinteractions(testtree(tree), trialtree(tree); kwargs...)`.
"""
function nearinteractions(tree, ::isBlockTree; kwargs...)
    return nearinteractions(testtree(tree), trialtree(tree); kwargs...)
end

"""
    nearinteractions(tree, trait; extractselfvalues=false, isnear=isnear)

Compute single-tree near interactions for the given tree trait.

When the tree is not balanced, ancestor-near leaves are included so near
interactions remain complete across different leaf levels.
"""
function nearinteractions(
    tree, ::A; extractselfvalues=false, isnear=isnear
) where {A<:AbstractTreeTrait}
    isleafnear = _LeafPredicateFunctor(isnear)

    selfv = Vector{Int}[]

    v = Vector{Int}[]
    nearvalues = Vector{Int}[]
    lk = Threads.SpinLock()

    isbalancedtree = checkbalancedtree(tree)

    @threads for node in leaves(tree)
        selfnearvalues = values(data(tree, node))

        lock(lk) do
            return push!(selfv, selfnearvalues)
        end

        nonselfnearvalues = Int[]

        isempty(selfnearvalues) && continue

        for nearnode in NearNodeIterator(tree, node; isnear=isnear)
            nearnode == node && continue
            appendvalues!(nonselfnearvalues, tree, nearnode)
        end

        if !isbalancedtree
            # for uniform trees, where all leaves are on same level we can skip this
            for parent in ParentUpwardsIterator(tree, node)
                for nearnode in NearNodeIterator(tree, parent; isnear=isleafnear)
                    appendvalues!(nonselfnearvalues, tree, nearnode)
                end
            end
        end

        isempty(nonselfnearvalues) && continue
        lock(lk) do
            push!(v, selfnearvalues)
            return push!(nearvalues, nonselfnearvalues)
        end
    end

    if extractselfvalues
        return selfv, v, nearvalues
    else
        prepend!(v, selfv)
        prepend!(nearvalues, selfv)
        return v, nearvalues
    end
end

"""
    nearinteractions(testtree, trialtree; extractselfvalues=false, isnear=isnear)

Compute two-tree near interactions.

Returns `(testvalues, trialvalues)`, where each `testvalues[i]` is paired with
`trialvalues[i]`. `extractselfvalues=true` is not supported for two-tree
interactions.
"""
function nearinteractions(testtree, trialtree; extractselfvalues=false, isnear=isnear)
    if extractselfvalues
        throw(
            ArgumentError(
                "extractselfvalues is only supported for single-tree nearinteractions"
            ),
        )
    end
    isleafnear = _LeafPredicateFunctor(isnear)

    testv = Vector{Int}[]
    trialv = Vector{Int}[]

    lk = Threads.SpinLock()
    arebalancedtrees = checkbalancedtree(testtree) && checkbalancedtree(trialtree)

    for testnode in leaves(testtree)
        testvalues = values(data(testtree, testnode))
        isempty(testvalues) && continue
        nearvalues = Int[]

        for nearnode in NearNodeIterator(trialtree, testtree, testnode; isnear=isnear)
            appendvalues!(nearvalues, trialtree, nearnode)
        end

        if !arebalancedtrees
            for parent in ParentUpwardsIterator(testtree, testnode)
                for nearnode in
                    NearNodeIterator(trialtree, testtree, parent; isnear=isleafnear)
                    appendvalues!(nearvalues, trialtree, nearnode)
                end
            end
        end

        isempty(nearvalues) && continue
        lock(lk) do
            push!(testv, testvalues)
            return push!(trialv, nearvalues)
        end
    end
    return testv, trialv
end
