
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

    elements, ad, _ = assemblydata(X; onlyactives=false)
    # @assert length(elements) == size(Σ, 2)

    gX = SimpleGraph(numfunctions(X))
    for e in Graphs.edges(gσ)
        elementa, elementb = src(e), dst(e)
        # J = length(ad[elementa])
        # I = length(ad[elementb])
        J = size(view(ad.data, :, :, elementa), 2)
        I = size(view(ad.data, :, :, elementb), 2)
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
    MetisForest(basis::BEAST.Space; graphweights=adjacencygraph(basis), builder=MetisForestBuilder())

Construct a `Forest` of `MetisTree`s from a `BEAST.Space`.

This method builds (or accepts) a function adjacency graph and graph weights for
`basis`, then delegates to the point/graph `MetisForest` constructor. The
builder carries the METIS construction settings, such as the number of requested
divisions per split.
"""
function MetisForest(
    basis::BEAST.Space;
    graphweights=adjacencygraph(basis),
    builder::MetisForestBuilder=MetisForestBuilder(),
)
    return MetisForest(positions(basis), graphweights[1], graphweights[2]; builder=builder)
end

"""
    buildforest(basis::BEAST.Space; graphweights=nothing, builder=MetisForestBuilder())

Build a METIS forest from a BEAST space using the builder workflow.

When `graphweights` is omitted, `adjacencygraph(basis)` supplies both the graph
and vertex weights. Passing explicit graph weights reuses them without
rebuilding the adjacency graph.
"""
function buildforest(
    basis::BEAST.Space;
    graphweights=nothing,
    builder::MetisForestBuilder=MetisForestBuilder(),
)
    graphweights = isnothing(graphweights) ? adjacencygraph(basis) : graphweights
    return MetisForest(basis; graphweights=graphweights, builder=builder)
end

"""
    MetisTree(basis::BEAST.Space; graphweights=adjacencygraph(basis), builder=MetisTreeBuilder())

Construct a `MetisTree` from a BEAST space using METIS graph partitioning.

This method builds (or accepts) a function adjacency graph and graph weights for
`basis`, then delegates to the point/graph `MetisTree` constructor. The builder
carries the METIS construction settings.
"""
function MetisTree(
    basis::BEAST.Space;
    graphweights=adjacencygraph(basis),
    builder::MetisTreeBuilder=MetisTreeBuilder(),
)
    return MetisTree(positions(basis), graphweights[1], graphweights[2]; builder=builder)
end

# Internal dispatch hook for the public
# buildtree(basis::BEAST.Space; builder=MetisTreeBuilder(), graphweights=nothing)
# route documented in H2BEASTTrees.jl.
function _buildbeasttree(basis::BEAST.Space, builder::MetisTreeBuilder, graphweights)
    graphweights = isnothing(graphweights) ? adjacencygraph(basis) : graphweights
    return MetisTree(basis; graphweights=graphweights, builder=builder)
end
