
function nearinteractions(tree; kwargs...)
    return nearinteractions(tree, treetrait(tree); kwargs...)
end

function nearinteractions(tree, ::isBlockTree; kwargs...)
    return nearinteractions(testtree(tree), trialtree(tree); kwargs...)
end

function nearinteractions(
    tree, treetrait::A; extractselfvalues=false, isnear=isnear
) where {A<:AbstractTreeTrait}
    isleafnear = _LeafNearFunctor(isnear)

    selfv = Vector{Int}[]

    v = Vector{Int}[]
    nearvalues = Vector{Int}[]
    lk = Threads.SpinLock()

    isbalancedtree = checkbalancedtree(tree)

    @threads for node in leaves(tree)
        selfnearvalues = values(tree, node)

        lock(lk) do
            return push!(selfv, selfnearvalues)
        end

        nonselfnearvalues = Int[]

        isempty(selfnearvalues) && continue

        for nearnode in NearNodeIterator(tree, node; isnear=isnear)
            nearnode == node && continue
            append!(nonselfnearvalues, values(tree, nearnode))
        end

        if !isbalancedtree
            # for uniform trees, where all leaves are on same level we can skip this
            for parent in ParentUpwardsIterator(tree, node)
                for nearnode in NearNodeIterator(tree, parent; isnear=isleafnear)
                    append!(nonselfnearvalues, values(tree, nearnode))
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

function nearinteractions(testtree, trialtree; extractselfvalues=false, isnear=isnear)
    # no selfvalues for two trees
    @assert !extractselfvalues
    isleafnear = _LeafNearFunctor(isnear)

    testv = Vector{Int}[]
    trialv = Vector{Int}[]

    lk = Threads.SpinLock()
    arebalancedtrees = checkbalancedtree(testtree) && checkbalancedtree(trialtree)

    for testnode in leaves(testtree)
        testvalues = values(testtree, testnode)
        isempty(testvalues) && continue
        nearvalues = Int[]

        for nearnode in NearNodeIterator(trialtree, testtree, testnode; isnear=isnear)
            append!(nearvalues, values(trialtree, nearnode))
        end

        if !arebalancedtrees
            for parent in ParentUpwardsIterator(testtree, testnode)
                for nearnode in
                    NearNodeIterator(trialtree, testtree, parent; isnear=isleafnear)
                    append!(nearvalues, values(trialtree, nearnode))
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
