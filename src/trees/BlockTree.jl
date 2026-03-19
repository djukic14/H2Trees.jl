
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
    testminvalues=0,
    trialminvalues=testminvalues,
    testmaxprotrusion=NaN,
    trialmaxprotrusion=testmaxprotrusion,
    testcomputeprotrusion=ComputeProtrusionFunctor(),
    trialcomputeprotrusion=ComputeProtrusionFunctor(),
)
    testcenter, testhalfsize = boundingbox(testpositions)
    trialcenter, trialhalfsize = boundingbox(trialpositions)

    minhalfsize, testroothalfsize, testminlevel, trialroothalfsize, trialminlevel, = adjusttwontreeblocktreeparameters(
        testhalfsize, trialhalfsize, minhalfsize
    )

    testtree = TwoNTree(
        SVector(testcenter...),
        testpositions,
        testroothalfsize,
        minhalfsize;
        minlevel=testminlevel,
        minvalues=testminvalues,
        maxprotrusion=testmaxprotrusion,
        computeprotrusion=testcomputeprotrusion,
    )

    trialtree = TwoNTree(
        SVector(trialcenter...),
        trialpositions,
        trialroothalfsize,
        minhalfsize;
        minlevel=trialminlevel,
        minvalues=trialminvalues,
        maxprotrusion=trialmaxprotrusion,
        computeprotrusion=trialcomputeprotrusion,
    )

    return BlockTree(testtree, trialtree)
end

function adjusttwontreeblocktreeparameters(testhalfsize, trialhalfsize, minhalfsize)
    if abs(testhalfsize - trialhalfsize) / max(testhalfsize, trialhalfsize) <
        eps(minhalfsize) * 1e6
        # case where both are the same size
        commonhalfsize = max(testhalfsize, trialhalfsize)

        return minhalfsize,
        roothalfsize(commonhalfsize, minhalfsize),
        1,
        roothalfsize(commonhalfsize, minhalfsize),
        1

    elseif trialhalfsize > testhalfsize
        # case where test tree is larger than trial tree
        testroothalfsize = roothalfsize(testhalfsize, minhalfsize)
        trialroothalfsize = roothalfsize(trialhalfsize, testroothalfsize)

        testminlevel = numberoflevels(trialroothalfsize, testroothalfsize) + 1
        trialminlevel = 1

        return minhalfsize, testroothalfsize, testminlevel, trialroothalfsize, trialminlevel
    else
        # case where trial tree is larger than test tree
        trialroothalfsize = roothalfsize(trialhalfsize, minhalfsize)
        testroothalfsize = roothalfsize(testhalfsize, trialroothalfsize)

        trialminlevel = numberoflevels(testroothalfsize, trialroothalfsize) + 1
        testminlevel = 1

        return minhalfsize, testroothalfsize, testminlevel, trialroothalfsize, trialminlevel
    end
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
