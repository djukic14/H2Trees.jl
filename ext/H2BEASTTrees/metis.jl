# function adjacencygraph(σ::BEAST.LagrangeBasies{0,-1,M,T,D,S}) where {M,T,D,S}
#     #TODO: generalize this notion of adjacency graph
#     #idea: two basis functions are connected in a graph if they have triangles in their
#     #support that touch each other
#     #implement: graph that tells us if two triangles touch

#     @assert all(x -> length(x) == 1, σ.fns)

#     edges = setdiff(skeleton(geometry(σ), 1), boundary(geometry(σ)))
#     Σ = connectivity(geometry(σ), edges, sign)
#     ΣΣ = Σ' * Σ

#     return Graph(ΣΣ)
# end

# two basis functions are connected in a graph if they have triangles in their support that touch each other
# this is the graph we want to partition with METIS
# weights for the graph are the areas of the triangles in the support of each basis function,
# divided by the number of basis functions supported on each triangle (so that the total weight of each triangle is equal to its area)

struct SignFunctor end

function (f::SignFunctor)(e)
    return sign(e)
end

function adjacencygraph(X)
    σ = lagrangecxd0(geometry(X))

    # we are abusing this as a graph that tells us if two triangles touch each other
    # we have to hope that BEAST internals do not change
    for (i, fn) in enumerate(σ.fns)
        @assert length(fn) == 1
        @assert fn[1].cellid == i
    end
    # edges = BEAST.setminus(BEAST.skeleton(geometry(σ), 1), BEAST.boundary(geometry(σ)))
    # Σ = BEAST.connectivity(geometry(σ), edges, SignFunctor())

    # this is not the actual Σ matrix because it does not work for open geometries,
    # but it works to get the adjacency graph (I hope so)
    edges = BEAST.skeleton(geometry(X), 1)
    faces = BEAST.skeleton(geometry(X), 2)
    Σ = BEAST.connectivity(faces, edges, SignFunctor())
    ΣΣ = Σ' * Σ
    gσ = Graph(ΣΣ) # graph that tells us if two triangles touch each other

    elements, ad, _ = assemblydata(X)

    @assert length(elements) == size(Σ, 2)

    gX = SimpleGraph(numfunctions(X))
    for e in Graphs.edges(gσ)
        elementa, elementb = src(e), dst(e)
        J = length(ad[elementa])
        I = length(ad[elementb])

        for j in 1:J
            for i in 1:I
                for (functionida, a) in ad[elementa][j]
                    iszero(a) && continue
                    for (functionidb, b) in ad[elementb][i]
                        iszero(b) && continue
                        functionida == functionidb && continue
                        add_edge!(gX, functionida, functionidb)
                    end
                end
            end
        end
    end

    graphweights = zeros(scalartype(X), numfunctions(X))

    functionsperelements = [length(ad[i]) for i in 1:length(elements)]

    for i in 1:numfunctions(X)
        for fn in X.fns[i]
            c = fn.cellid
            graphweights[i] +=
                BEAST.volume(BEAST.chart(geometry(X), c)) / functionsperelements[c]
        end
    end
    return gX, graphweights
end

function MetisForest(
    basis::BEAST.Space,
    numdivisions::Int;
    splitunconnectedpartitions=false,
    minlevel::Int=1,
    minvalues::Int=numdivisions,
    root::Int=1,
    balldata=BoundingBallData,
    graphweights=adjacencygraph(basis),
)
    return MetisForest(
        positions(basis),
        graphweights[1],
        graphweights[2],
        numdivisions;
        splitunconnectedpartitions=splitunconnectedpartitions,
        minlevel=minlevel,
        minvalues=minvalues,
        root=root,
        balldata=balldata,
    )
end

function MetisTree(
    basis::BEAST.Space,
    numdivisions::Int;
    splitunconnectedpartitions=false,
    minlevel::Int=1,
    minvalues::Int=numdivisions,
    root::Int=1,
    balldata=BoundingBallData,
    graphweights=adjacencygraph(basis),
)
    return MetisTree(
        positions(basis),
        graphweights[1],
        graphweights[2],
        numdivisions;
        splitunconnectedpartitions=splitunconnectedpartitions,
        minlevel=minlevel,
        minvalues=minvalues,
        root=root,
        balldata=balldata,
    )
end
