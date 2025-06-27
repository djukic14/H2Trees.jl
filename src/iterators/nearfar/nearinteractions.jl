
function nearinteractions(tree; kwargs...)
    return nearinteractions(tree, treetrait(tree); kwargs...)
end

function nearinteractions(tree, ::isBlockTree; kwargs...)
    !arepointsunique(tree) &&
        return nearinteractionsnonunique(testtree(tree), trialtree(tree); kwargs...)
    return nearinteractions(testtree(tree), trialtree(tree); kwargs...)
end

function nearinteractions(
    tree, treetrait::A; extractselfvalues=false, isnear=isnear
) where {A<:AbstractTreeTrait}
    !arepointsunique(tree) && return nearinteractionsnonunique(
        tree, treetrait; extractselfvalues=extractselfvalues, isnear=isnear
    )

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

function setappend!(set, values)
    for val in values
        push!(set, val)
    end
end

function nearinteractionsnonunique(
    tree, treetrait::A; extractselfvalues=false, isnear=isnear
) where {A<:AbstractTreeTrait}
    valatnodes = valuesatnodes(tree)
    natvalues = nodesatvalues(tree, valatnodes)

    selfv = Vector{Int}[]

    v = Vector{Int}[]
    nearvalues = Vector{Int}[]
    lk = Threads.SpinLock()

    @assert checkbalancedtree(tree)

    @threads for node in leaves(tree)
        selfnearvalues = Int[]
        for val in values(tree, node)
            # values that are in multiple boxes are treated later
            length(valatnodes[val]) > 1 && continue
            push!(selfnearvalues, val)
        end

        isempty(selfnearvalues) && continue

        lock(lk) do
            return push!(selfv, selfnearvalues)
        end

        nonselfnearvalues = Set{Int}() # make sure that nearvalues are unique

        for nearnode in NearNodeIterator(tree, node; isnear=isnear)
            nearnode == node && continue
            setappend!(nonselfnearvalues, values(tree, nearnode))
        end

        nonselfnearvalues = collect(nonselfnearvalues)
        isempty(nonselfnearvalues) && continue
        lock(lk) do
            push!(v, selfnearvalues)
            return push!(nearvalues, nonselfnearvalues)
        end
    end

    @threads for (boxes, sharedvalues) in collect(natvalues)
        ((length(boxes) <= 1) || (isempty(sharedvalues))) && continue
        # consider values that are in multiple boxes

        @lock lk push!(selfv, sharedvalues)

        _nearvalues = Set{Int}()
        for box in boxes
            for node in H2Trees.NearNodeIterator(tree, box; isnear=isnear)
                for val in H2Trees.values(tree, node)
                    val in sharedvalues && continue
                    push!(_nearvalues, val)
                end
            end
        end

        _nearvalues = collect(_nearvalues)
        isempty(_nearvalues) && continue

        lock(lk) do
            push!(v, sharedvalues)
            return push!(nearvalues, _nearvalues)
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

function nearinteractionsnonunique(
    testtree, trialtree; extractselfvalues=false, isnear=isnear
)
    # no selfvalues for two trees
    @assert !extractselfvalues

    valatnodes = valuesatnodes(testtree)
    natvalues = nodesatvalues(testtree, valatnodes)

    testv = Vector{Int}[]
    trialv = Vector{Int}[]

    lk = Threads.SpinLock()
    @assert checkbalancedtree(testtree) && checkbalancedtree(trialtree)

    @threads for testnode in leaves(testtree)
        testvalues = Int[]
        for val in values(testtree, testnode)
            # values that are in multiple boxes are treated later
            length(valatnodes[val]) > 1 && continue
            push!(testvalues, val)
        end

        isempty(testvalues) && continue

        nearvalues = Set{Int}()

        for nearnode in NearNodeIterator(trialtree, testtree, testnode; isnear=isnear)
            setappend!(nearvalues, values(trialtree, nearnode))
        end

        nearvalues = collect(nearvalues)
        isempty(nearvalues) && continue
        lock(lk) do
            push!(testv, testvalues)
            return push!(trialv, nearvalues)
        end
    end

    @threads for (boxes, sharedvalues) in collect(natvalues)
        ((length(boxes) <= 1) || (isempty(sharedvalues))) && continue

        nearvalues = Set{Int}()
        for box in boxes
            for nearnode in
                H2Trees.NearNodeIterator(trialtree, testtree, box; isnear=isnear)
                for val in H2Trees.values(trialtree, nearnode)
                    push!(nearvalues, val)
                end
            end
        end
        nearvalues = collect(nearvalues)
        isempty(nearvalues) && continue

        lock(lk) do
            push!(testv, sharedvalues)
            return push!(trialv, nearvalues)
        end
    end
    return testv, trialv
end
