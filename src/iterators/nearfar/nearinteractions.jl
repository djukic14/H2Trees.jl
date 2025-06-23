
function nearinteractions(tree; kwargs...)
    return nearinteractions(tree, treetrait(tree); kwargs...)
end

function nearinteractions(tree, ::isBlockTree; kwargs...)
    !uniquepointstree(tree) &&
        return nearinteractionsnonunique(testtree(tree), trialtree(tree); kwargs...)
    return nearinteractions(testtree(tree), trialtree(tree); kwargs...)
end

function nearinteractions(
    tree,
    treetrait::A;
    size=typemax(Int),
    minsize=nothing,
    extractselfvalues=false,
    isnear=isnear,
) where {A<:AbstractTreeTrait}
    !uniquepointstree(tree) && return nearinteractionsnonunique(
        tree,
        treetrait;
        size=size,
        minsize=minsize,
        extractselfvalues=extractselfvalues,
        isnear=isnear,
    )

    isleafnear = LeafNearFunctor(isnear)

    selfv = Vector{Int}[]

    v = Vector{Int}[]
    nearvalues = Vector{Int}[]
    lk = Threads.SpinLock()

    isbalancedtree = checkbalancedtree(tree)

    @threads for node in leaves(tree)
        selfnearvalues = values(tree, node)

        lock(lk) do
            for (i, is) in enumerate(chunks(selfnearvalues; size=size, minsize=minsize))
                isempty(is) && continue
                for (j, js) in enumerate(chunks(selfnearvalues; size=size, minsize=minsize))
                    isempty(js) && continue
                    if i == j
                        push!(selfv, is)
                    else
                        push!(v, is)
                        push!(nearvalues, js)
                    end
                end
            end
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
            for selfnearvalues in chunks(selfnearvalues; size=size, minsize=minsize)
                for nonselfnearvalues in
                    chunks(nonselfnearvalues; size=size, minsize=minsize)
                    push!(v, selfnearvalues)
                    push!(nearvalues, nonselfnearvalues)
                end
            end
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

function nearinteractions(
    testtree, trialtree; size=500, minsize=nothing, extractselfvalues=false, isnear=isnear
)
    # no selfvalues for two trees
    @assert !extractselfvalues
    isleafnear = LeafNearFunctor(isnear)

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
            for testvalues in chunks(testvalues; size=size, minsize=minsize)
                for nearvalues in chunks(nearvalues; size=size, minsize=minsize)
                    push!(testv, testvalues)
                    push!(trialv, nearvalues)
                end
            end
        end
    end
    return testv, trialv
end

function isprimarybox(nodes, node)
    return nodes[begin] == node
end

function setappend!(set, values)
    for val in values
        push!(set, val)
    end
end

function nearinteractionsnonunique(
    tree,
    treetrait::A;
    size=typemax(Int),
    minsize=nothing,
    extractselfvalues=false,
    isnear=isnear,
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
            for (i, is) in enumerate(chunks(selfnearvalues; size=size, minsize=minsize))
                isempty(is) && continue
                for (j, js) in enumerate(chunks(selfnearvalues; size=size, minsize=minsize))
                    isempty(js) && continue
                    if i == j
                        push!(selfv, is)
                    else
                        push!(v, is)
                        push!(nearvalues, js)
                    end
                end
            end
        end

        nonselfnearvalues = Set{Int}() # make sure that nearvalues are unique

        for nearnode in NearNodeIterator(tree, node; isnear=isnear)
            nearnode == node && continue
            setappend!(nonselfnearvalues, values(tree, nearnode))
        end

        nonselfnearvalues = collect(nonselfnearvalues)
        isempty(nonselfnearvalues) && continue
        lock(lk) do
            for selfnearvalues in chunks(selfnearvalues; size=size, minsize=minsize)
                for nonselfnearvalues in
                    chunks(nonselfnearvalues; size=size, minsize=minsize)
                    push!(v, selfnearvalues)
                    push!(nearvalues, nonselfnearvalues)
                end
            end
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
            for sharedvalues in chunks(sharedvalues; size=size, minsize=minsize)
                for _nearvalues in chunks(_nearvalues; size=size, minsize=minsize)
                    push!(v, sharedvalues)
                    push!(nearvalues, _nearvalues)
                end
            end
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
    testtree, trialtree; size=500, minsize=nothing, extractselfvalues=false, isnear=isnear
)
    @warn "Bug here"
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
            for testvalues in chunks(testvalues; size=size, minsize=minsize)
                for nearvalues in chunks(nearvalues; size=size, minsize=minsize)
                    push!(testv, testvalues)
                    push!(trialv, nearvalues)
                end
            end
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
            for sharedvalues in chunks(sharedvalues; size=size, minsize=minsize)
                for nearvalues in chunks(nearvalues; size=size, minsize=minsize)
                    push!(testv, sharedvalues)
                    push!(trialv, nearvalues)
                end
            end
        end
    end
    return testv, trialv
end
