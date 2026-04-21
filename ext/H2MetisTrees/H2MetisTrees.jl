module H2MetisTrees

using StaticArrays
using Graphs
using Metis
using Metis.LibMetis: idx_t, @check, METIS_PartGraphKway, METIS_PartGraphRecursive

import H2Trees: metispartition

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

"""
        partition(G::Metis.Graph, nparts::Integer; alg=:KWAY, options=metisoptions, vertexweights=true)

Partition a METIS graph into `nparts` parts.

This is a thin wrapper around METIS partitioning routines that supports both
recursive bisection (`:RECURSIVE`) and k-way partitioning (`:KWAY`).

# Arguments

    - `G::Metis.Graph`: Input graph in METIS format.
    - `nparts::Integer`: Number of requested partitions.
    - `alg`: Partitioning algorithm (`:KWAY` or `:RECURSIVE`, default: `:KWAY`).
    - `options`: METIS options vector (default: `metisoptions`).
    - `vertexweights`: Whether to use `G.vwgt` as vertex weights (default: `true`).

# Returns

A vector of partition labels with one entry per vertex.
"""
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

"""
        computemetisweights(::Type{T}, weights, targetmax=1000) where {T}

Convert floating-point-like weights into positive integer METIS vertex weights.

The input is transformed with `log1p`, normalized to `[0, 1]`, scaled by
`targetmax`, and rounded to integers. The output is clamped to be at least `1`,
as required by METIS.

# Arguments

    - `::Type{T}`: Target integer type for METIS weights.
    - `weights`: Input vertex weights.
    - `targetmax`: Maximum scale used before rounding (default: `1000`).

# Returns

A vector of type `T` containing strictly positive integer weights.
"""
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

"""
        metispartition(G::Graphs.SimpleGraph, vertexweights, numberofdivisions::Int; splitunconnectedpartitions::Bool=false, alg=:KWAY)

Partition a simple graph using METIS with optional post-processing of disconnected parts.

The function converts `G` to a METIS graph, computes integer vertex weights,
and runs weighted partitioning first. If METIS returns a single part, it retries
without vertex weights using `fallbackmetisoptions`. Optionally, each resulting
part can be split into connected components.

# Arguments

    - `G::Graphs.SimpleGraph`: Input graph to partition.
    - `vertexweights`: Vertex weights used to bias partitioning.
    - `numberofdivisions::Int`: Number of requested partitions.
    - `splitunconnectedpartitions::Bool`: Split each METIS part into connected components (default: `false`).
    - `alg`: Partitioning algorithm to use (`:KWAY` or `:RECURSIVE`, default: `:KWAY`).

# Returns

A vector of partition labels with one entry per vertex.
"""
function metispartition(
    G::Graphs.SimpleGraph,
    vertexweights,
    numberofdivisions::Int;
    splitunconnectedpartitions::Bool=false,
    alg=:KWAY,
)
    metisG = Metis.graph(G)
    metisG = Metis.Graph(
        metisG.nvtxs,
        metisG.xadj,
        metisG.adjncy,
        computemetisweights(typeof(metisG.nvtxs), vertexweights),
        metisG.adjwgt,
    )

    part = partition(
        metisG, numberofdivisions; options=metisoptions, vertexweights=true, alg=alg
    )

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
        newpart = zeros(Int, length(part))
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

end # module H2MetisTrees
