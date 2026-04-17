
struct SignFunctor end

function (f::SignFunctor)(e)
    return sign(e)
end

"""
    adjacencygraph(X::BEAST.Space)

Construct an adjacency graph and corresponding vertex weights for a BEAST space.

Two basis functions are considered adjacent if they are supported on touching mesh
elements. The graph weights are assembled from element areas, distributed across
all basis functions supported on each element.

# Arguments

  - `X::BEAST.Space`: Basis function space used to build adjacency and weights.

# Returns

A tuple `(graph, graphweights)` where:

  - `graph::SimpleGraph` is the adjacency graph between basis functions.
  - `graphweights` contains one weight per basis function.

# See also

`MetisTree`, `MetisForest`.
"""
function adjacencygraph(X::BEAST.Space)
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

    # @assert length(elements) == size(Σ, 2)

    gX = SimpleGraph(numfunctions(X))
    for e in Graphs.edges(gσ)
        elementa, elementb = src(e), dst(e)
        J = length(ad[elementa])
        I = length(ad[elementb])
        for j in 1:J
            for i in 1:I
                for (functionida, a) in view(ad.data[:, j, elementa])
                    iszero(a) && continue
                    for (functionidb, b) in view(ad.data[:, i, elementb])
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

"""
    MetisForest(basis::BEAST.Space, numdivisions::Int; graphweights=adjacencygraph(basis), kwargs...)

Construct a `Forest` of `MetisTree`s from a `BEAST.Space`.

This method builds (or accepts) a function adjacency graph and graph weights for
`basis`, then forwards to the point/graph `MetisForest` constructor.

# Arguments

  - `basis::BEAST.Space`: Function space used to build the partitioning graph.
  - `numdivisions::Int`: Number of partitions requested at each split.
  - `graphweights`: Tuple `(graph, weights)` used for partitioning (default: `adjacencygraph(basis)`).
  - `kwargs...`: Additional keyword arguments forwarded to `MetisForest(points, graph, weights, ...)`.

# Returns

A `Forest` containing one `MetisTree` per connected component of the basis graph.
"""
function MetisForest(
    basis::BEAST.Space, numdivisions::Int; graphweights=adjacencygraph(basis), kwargs...
)
    return MetisForest(
        positions(basis), graphweights[1], graphweights[2], numdivisions; kwargs...
    )
end

"""
    MetisTree(basis::BEAST.Space, numdivisions::Int; graphweights=adjacencygraph(basis), kwargs...)

Construct a `BoundingBallTree` from a `BEAST.Space` using METIS graph partitioning.

This method builds (or accepts) a function adjacency graph and graph weights for
`basis`, then forwards to the point/graph `MetisTree` constructor.

# Arguments

  - `basis::BEAST.Space`: Function space used to build the partitioning graph.
  - `numdivisions::Int`: Number of partitions requested at each split.
  - `graphweights`: Tuple `(graph, weights)` used for partitioning (default: `adjacencygraph(basis)`).
  - `kwargs...`: Additional keyword arguments forwarded to `MetisTree(points, graph, weights, ...)`.

# Returns

A `BoundingBallTree` with basis-function positions organized hierarchically.
"""
function MetisTree(
    basis::BEAST.Space, numdivisions::Int; graphweights=adjacencygraph(basis), kwargs...
)
    return MetisTree(
        positions(basis), graphweights[1], graphweights[2], numdivisions; kwargs...
    )
end
