
struct BlockTree{T}
    testcluster::T
    trialcluster::T
end

function TwoNTree(
    testpositions,
    trialpositions,
    minhalfsize;
    safetyfactorboundingbox=1,
    pointsunique=H2Trees.UniquePoints(),
)
    testcenter, testhalfsize = boundingbox(testpositions)
    trialcenter, trialhalfsize = boundingbox(trialpositions)

    trialnlevels = numberoflevels(trialhalfsize * safetyfactorboundingbox, minhalfsize)
    testnlevels = numberoflevels(testhalfsize * safetyfactorboundingbox, minhalfsize)

    minleveltrial = trialnlevels >= testnlevels ? 1 : testnlevels - trialnlevels + 1
    minleveltest = trialnlevels >= testnlevels ? trialnlevels - testnlevels + 1 : 1

    testtree = TwoNTree(
        SVector(testcenter...),
        testpositions,
        roothalfsize(minhalfsize, testnlevels),
        minhalfsize;
        minlevel=minleveltest,
        pointsunique=pointsunique,
    )

    trialtree = TwoNTree(
        SVector(trialcenter...),
        trialpositions,
        roothalfsize(minhalfsize, trialnlevels),
        minhalfsize;
        minlevel=minleveltrial,
        pointsunique=pointsunique,
    )
    #TODO: create both trees with level 1 and then update level afterwards --> important for minvalues ≠ 0

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
    return Base.eltype(H2Trees.testtree(tree))
end
