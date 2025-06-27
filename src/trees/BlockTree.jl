
"""
    BlockTree
"""
struct BlockTree{T}
    testcluster::T
    trialcluster::T
end

function TwoNTree(
    testpositions,
    trialpositions,
    minhalfsize;
    safetyfactorboundingbox=1,
    minvaluestest=0,
    minvaluestrial=0,
)
    # We cannot just create the trees from the minimum number of values, because
    # the boxes at the same level across both trees have to be the same size.
    @assert minhalfsize > 0 "Minimum halfsize must be greater than zero."

    testcenter, testhalfsize = boundingbox(testpositions)
    trialcenter, trialhalfsize = boundingbox(trialpositions)

    trialnlevels = numberoflevels(trialhalfsize * safetyfactorboundingbox, minhalfsize)
    testnlevels = numberoflevels(testhalfsize * safetyfactorboundingbox, minhalfsize)

    minleveltrial = trialnlevels >= testnlevels ? 1 : testnlevels - trialnlevels + 1
    minleveltest = trialnlevels >= testnlevels ? trialnlevels - testnlevels + 1 : 1

    testtree = TwoNTree(
        SVector(testcenter...),
        testpositions,
        roothalfsize(testhalfsize, minhalfsize),
        minhalfsize;
        minlevel=minleveltest,
        minvalues=minvaluestest,
    )

    trialtree = TwoNTree(
        SVector(trialcenter...),
        trialpositions,
        roothalfsize(trialhalfsize, minhalfsize),
        minhalfsize;
        minlevel=minleveltrial,
        minvalues=minvaluestrial,
    )

    return BlockTree(testtree, trialtree)
end

function treetrait(::Type{BlockTree{T}}) where {T}
    return isBlockTree()
end

function testtree(tree::BlockTree) # if relevant
    return tree.testcluster
end

function trialtree(tree::BlockTree) # if relevant
    return tree.trialcluster
end

function Base.eltype(tree::BlockTree)
    return Base.eltype(testtree(tree))
end
