function QuadPointsTree(
    space::B,
    minhalfsize;
    adata=assemblydata(space),
    qstrat=BEAST.SingleNumQStrat(5),
    qdata=BEAST.quadpoints(refspace(space), space, (qstrat.quad_rule,)),
    positions=nothing,
) where {B<:BEASTnurbs.BsplineBasis}
    positions = if isnothing(positions)
        pos = SVector{3,Float64}[] #TODO: make this typestable
        for i in eachindex(qdata)
            for q in qdata[i]
                push!(pos, BEAST.cartesian(q.point))
            end
        end
        convert(Vector{typeof(pos[1])}, pos)
    end

    tree = TwoNTree(positions, minhalfsize; pointsunique=H2Trees.NonUniquePoints())

    cU, cV = BEASTnurbs.numBezierCells(space)

    qpointsindices = CartesianIndices((length(first(qdata)), size(qdata)...))

    for leaf in H2Trees.leaves(tree)
        vals = H2Trees.values(tree, leaf)
        newvals = Set{Int}()

        for val in vals
            _, pInd, cellU, cellV = qpointsindices[val].I
            cell = NURBS.cellCart2Lin(cellU, cellV, pInd, cU, cV)

            for (functionid, _) in adata[cell]
                push!(newvals, functionid)
            end
        end

        empty!(H2Trees.values(tree, leaf))
        append!(H2Trees.values(tree, leaf), sort!(collect(newvals)))
    end

    return tree
end
