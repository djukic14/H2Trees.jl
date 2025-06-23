function _quadpointpositions(qdata)
    positions = []
    for i in eachindex(qdata)
        for q in qdata[i]
            push!(positions, BEAST.cartesian(q.point))
        end
    end
    positions = convert(Vector{typeof(positions[1])}, positions)
    return positions
end

function QuadPointsTree(
    space::BEAST.Space,
    minhalfsize;
    adinfo=assemblydata(space),
    qstrat=BEAST.SingleNumQStrat(5),
    qdata=BEAST.quadpoints(refspace(space), adinfo[1], (qstrat.quad_rule,)),
    positions=_quadpointpositions(qdata),
    minlevel=1,
)
    els, adata, _ = adinfo
    ref = refspace(space)

    tree = TwoNTree(
        positions, minhalfsize; pointsunique=NonUniquePoints(), minlevel=minlevel
    )
    qpointsindices = CartesianIndices((length(first(qdata)), length(els)))

    for leaf in H2Trees.leaves(tree)
        vals = H2Trees.values(tree, leaf)
        newvals = Set{Int}()
        for i in 1:numfunctions(ref)
            for val in vals
                for (functionid, _) in adata[qpointsindices[val].I[2], i]
                    push!(newvals, functionid)
                end
            end
        end

        empty!(H2Trees.values(tree, leaf))
        append!(H2Trees.values(tree, leaf), sort!(collect(newvals)))
    end

    return tree
end

function QuadPointsTree(
    testspace::BEAST.Space,
    trialspace::BEAST.Space,
    minhalfsize;
    testadinfo=assemblydata(testspace),
    trialadinfo=assemblydata(trialspace),
    testqstrat=BEAST.SingleNumQStrat(5),
    trialqstrat=BEAST.SingleNumQStrat(5),
    testqdata=BEAST.quadpoints(refspace(testspace), testadinfo[1], (testqstrat.quad_rule,)),
    trialqdata=BEAST.quadpoints(
        refspace(trialspace), trialadinfo[1], (trialqstrat.quad_rule,)
    ),
)
    testpositions = _quadpointpositions(testqdata)
    trialpositions = _quadpointpositions(trialqdata)

    #TODO: remove code duplication from BlockTree
    _, testhalfsize = boundingbox(testpositions)
    _, trialhalfsize = boundingbox(trialpositions)

    trialnlevels = numberoflevels(trialhalfsize, minhalfsize)
    testnlevels = numberoflevels(testhalfsize, minhalfsize)

    minleveltrial = trialnlevels >= testnlevels ? 1 : testnlevels - trialnlevels + 1
    minleveltest = trialnlevels >= testnlevels ? trialnlevels - testnlevels + 1 : 1

    testtree = QuadPointsTree(
        testspace,
        minhalfsize;
        adinfo=testadinfo,
        qstrat=testqstrat,
        qdata=testqdata,
        positions=testpositions,
        minlevel=minleveltest,
    )

    trialtree = QuadPointsTree(
        trialspace,
        minhalfsize;
        adinfo=trialadinfo,
        qstrat=trialqstrat,
        qdata=trialqdata,
        positions=trialpositions,
        minlevel=minleveltrial,
    )

    return BlockTree(testtree, trialtree)
end
