module H2MetisTrees

using BoundingSphere
using StaticArrays
using Graphs
using Metis
using Metis.LibMetis: idx_t, @check, METIS_PartGraphKway
using H2Trees
import H2Trees:
    BoundingBallTree,
    BoundingBallData,
    values,
    data,
    lastnode,
    level,
    _adjustnodesatlevels!,
    updateradii!,
    _updatechild!,
    _updatenextsibling!
"""
    fallbackmetisoptions

METIS options vector with minimal configuration.

This option set enables:

  - `METIS_OPTION_CONTIG = 1`: forces each partition to be contiguous
    (i.e., each block forms a connected subgraph).

All other METIS options remain unset (`-1`), meaning METIS will use defaults.
"""
const fallbackmetisoptions = begin
    _options = fill(Cint(-1), Metis.METIS_NOPTIONS)
    _options[Int(Metis.METIS_OPTION_CONTIG) + 1] = 1
    SVector{Metis.METIS_NOPTIONS}(_options)
end

"""
    metisoptions

METIS options vector with extended configuration.

This option set enables:

  - `METIS_OPTION_CONTIG = 1`: ensures partitions are contiguous
  - `METIS_OPTION_NUMBERING = 1`: switches METIS to 1-based numbering

All other options remain at default (`-1`).
"""
const metisoptions = begin
    _options = fill(Cint(-1), Metis.METIS_NOPTIONS)
    _options[Int(Metis.METIS_OPTION_CONTIG) + 1] = 1
    _options[Int(Metis.METIS_OPTION_NUMBERING) + 1] = 1
    SVector{Metis.METIS_NOPTIONS}(_options)
end

function partition(
    G::Metis.Graph, nparts::Integer; alg=:KWAY, options=metisoptions, vertexweights=true
)
    edgecut = fill(idx_t(0), 1)

    part = Vector{idx_t}(undef, G.nvtxs)

    weights = vertexweights ? G.vwgt : C_NULL

    nparts == 1 && return fill!(part, 1) # https://github.com/JuliaSparse/Metis.jl/issues/49
    if alg === :RECURSIVE
        @check METIS_PartGraphRecursive(
            Ref{idx_t}(G.nvtxs),
            Ref{idx_t}(1),
            G.xadj,
            G.adjncy,
            weights,
            C_NULL,
            G.adjwgt,
            Ref{idx_t}(nparts),
            C_NULL,
            C_NULL,
            options,
            edgecut,
            part,
        )
    elseif alg === :KWAY
        @check METIS_PartGraphKway(
            Ref{idx_t}(G.nvtxs),
            Ref{idx_t}(1),
            G.xadj,
            G.adjncy,
            weights,
            C_NULL,
            G.adjwgt,
            Ref{idx_t}(nparts),
            C_NULL,
            C_NULL,
            options,
            edgecut,
            part,
        )
    else
        throw(ArgumentError("unknown algorithm $(repr(alg))"))
    end
    return part
end

function computemetisweights(::Type{T}, weights, targetmax=1000) where {T}
    w = Float64.(weights)

    if all(iszero, w)
        return fill(T(1), length(w))
    end
    w = Float64.(weights)
    w = log1p.(w)
    w ./= maximum(w) # normalize to [0, 1]
    vwgt = max.(1, round.(Int, w .* targetmax))

    # ensure no zeros
    vwgt[vwgt .< 1] .= 1

    return T.(vwgt)
end

function metispartition(
    G::Graphs.SimpleGraph,
    vertexweights,
    numberofdivisions::Int;
    splitunconnectedpartitions::Bool=false,
)
    metisG = Metis.graph(G)
    metisG = Metis.Graph(
        metisG.nvtxs,
        metisG.xadj,
        metisG.adjncy,
        computemetisweights(typeof(metisG.nvtxs), vertexweights),
        metisG.adjwgt,
    )

    part = partition(metisG, numberofdivisions; options=metisoptions, vertexweights=true)

    if maximum(part) == minimum(part)
        @warn "METIS partitioning failed to produce multiple parts, retrying without vertex weights"
        part = partition(
            metisG, numberofdivisions; options=fallbackmetisoptions, vertexweights=false
        )
    end

    if splitunconnectedpartitions
        @warn "Splitting unconnected partitions into connected components."
        partitions = [Vector{Int}() for _ in 1:maximum(part)]
        for (i, p) in enumerate(part)
            push!(partitions[p], i)
        end

        partitionid = 0
        newpart = zeros(Vector{Int}, 0)
        for p in eachindex(partitions)
            subgraph, localtoglobal = induced_subgraph(G, partitions[p])
            components = connected_components(subgraph)

            for (i, comp) in enumerate(components)
                newpart[localtoglobal[comp]] .= partitionid + i
            end
            partitionid += length(components)
        end
        part = newpart
    end

    uniqueparts = unique(part)
    if length(uniqueparts) > numberofdivisions # TODO: is this really an issue or should we just warn?
        @warn "METIS partitioning produced too many partitions, got $(length(uniqueparts)) partitions instead of $numberofdivisions"
    end

    return part
end

function MetisTree(
    points::AbstractVector{SVector{N,T}},
    graph,
    pointgraphweights,
    numdivisions::Int;
    splitunconnectedpartitions=false,
    minlevel::Int=1,
    minvalues::Int=numdivisions,
    root::Int=1,
    balldata=BoundingBallData,
) where {N,T}
    center, radius = boundingsphere(points)

    #TODO: return vector of trees: one for each connected component of the graph
    tree = BoundingBallTree(center, radius; minlevel=minlevel, root=root, balldata=balldata)
    append!(values(data(tree, root)), collect(1:length(points)))

    splitnode!(
        tree,
        points,
        pointgraphweights,
        graph,
        root,
        numdivisions;
        splitunconnectedpartitions=splitunconnectedpartitions,
        minvalues=minvalues,
        balldata=balldata,
    )
    _adjustnodesatlevels!(tree)
    # updateradii!(tree; update=updateradii)
    return tree
end

function splitnode!(
    tree::BoundingBallTree,
    points::AbstractVector{SVector{N,T}},
    pointgraphweights,
    graph,
    node,
    numdivisions::Int;
    minvalues::Int=numdivisions,
    splitunconnectedpartitions=false,
    balldata=BoundingBallData,
) where {N,T}
    length(values(tree, node)) <= max(minvalues, numdivisions) && return tree

    subgraph, _ = induced_subgraph(graph, values(tree, node))
    part = metispartition(subgraph, pointgraphweights[values(tree, node)], numdivisions)

    partitions = [Vector{Int}() for _ in 1:maximum(part)]
    for (i, part) in enumerate(part)
        push!(partitions[part], values(tree, node)[i])
    end

    centers = Vector{SVector{N,T}}(undef, length(partitions))
    radii = Vector{T}(undef, length(partitions))

    for i in eachindex(partitions)
        centers[i], radii[i] = boundingsphere(points[partitions[i]])
    end

    _updatechild!(tree, node, lastnode(tree) + 1)

    for i in eachindex(partitions)
        dat = balldata(partitions[i], centers[i], radii[i], level(tree, node) + 1)
        childnode = lastnode(tree) + 1
        push!(tree.nodes, Node(dat, 0, node, 0))
        splitnode!(
            tree,
            points,
            pointgraphweights,
            graph,
            childnode,
            numdivisions;
            minvalues=minvalues,
            splitunconnectedpartitions=splitunconnectedpartitions,
            balldata=balldata,
        )
        _updatenextsibling!(
            tree, childnode, i == last(eachindex(centers)) ? 0 : lastnode(tree) + 1
        )
    end
    empty!(values(data(tree, node)))

    return tree
end

end # module H2MetisTrees
