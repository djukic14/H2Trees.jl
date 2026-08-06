module H2MetisTrees

using StaticArrays
using Graphs
using Metis
using Metis.LibMetis:
    idx_t,
    @check,
    METIS_PartGraphKway,
    METIS_PartGraphRecursive,
    METIS_OPTION_SEED,
    MetisError

import H2Trees: metispartition

# The same seed `KMeansTreeBuilder` uses by default (see its docstring): both partitioning
# strategies need a fixed seed so results are reproducible run to run and platform to platform,
# rather than depending on METIS's own default (time-based on at least some builds) seeding.
# Without this, `METIS_OPTION_CONTIG = 1` below can fail non-deterministically -- METIS refuses to
# produce a contiguous partition once its (seed-dependent) recursive splitting reaches a subgraph
# that is already disconnected, and an unset seed makes that reachable on some platforms/builds
# and not others for the exact same input graph.
const metisseed = 1234

"""
    fallbackmetisoptions

METIS options used when weighted partitioning fails to split the graph.

The fallback keeps `METIS_OPTION_CONTIG = 1` so each partition forms a
connected subgraph, and fixes `METIS_OPTION_SEED` for reproducibility. All
other options remain unset, letting METIS choose its defaults.
"""
const fallbackmetisoptions = begin
    _options = fill(Cint(-1), Metis.METIS_NOPTIONS)
    _options[Int(Metis.METIS_OPTION_CONTIG) + 1] = 1
    _options[Int(METIS_OPTION_SEED) + 1] = metisseed
    SVector{Metis.METIS_NOPTIONS}(_options)
end

"""
    metisoptions

Default METIS options used by H2Trees.

The options request contiguous partitions, Julia-style 1-based numbering,
and a fixed `METIS_OPTION_SEED` for reproducibility. All other options
remain unset, letting METIS choose its defaults.
"""
const metisoptions = begin
    _options = fill(Cint(-1), Metis.METIS_NOPTIONS)
    _options[Int(Metis.METIS_OPTION_CONTIG) + 1] = 1
    _options[Int(Metis.METIS_OPTION_NUMBERING) + 1] = 1
    _options[Int(METIS_OPTION_SEED) + 1] = metisseed
    SVector{Metis.METIS_NOPTIONS}(_options)
end

"""
    relaxedmetisoptions

METIS options used when a contiguous partition is infeasible.

`METIS_OPTION_CONTIG` requiring every partition to be internally connected
also requires the *input* graph itself to already be connected; on a
disconnected input graph (e.g. a mesh made of several physically separate
pieces), some METIS builds return an error for this instead of silently
falling back, so [`metispartition`](@ref) drops the contiguity constraint
entirely as a last resort. `METIS_OPTION_SEED` is still fixed for
reproducibility.
"""
const relaxedmetisoptions = begin
    _options = fill(Cint(-1), Metis.METIS_NOPTIONS)
    _options[Int(Metis.METIS_OPTION_NUMBERING) + 1] = 1
    _options[Int(METIS_OPTION_SEED) + 1] = metisseed
    SVector{Metis.METIS_NOPTIONS}(_options)
end

"""
    partition(G::Metis.Graph, nparts::Integer; alg=:KWAY, options=metisoptions, vertexweights=true)

Partition a METIS graph into `nparts` parts.

This thin wrapper dispatches to either k-way partitioning (`alg=:KWAY`) or
recursive bisection (`alg=:RECURSIVE`). Set `vertexweights=false` to ignore
`G.vwgt`. The returned vector contains one partition label per vertex.
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
    metispartition(G::Graphs.SimpleGraph, vertexweights, numberofdivisions::Int; splitunconnectedpartitions=false, alg=:KWAY)

Partition a simple graph using METIS with optional post-processing of disconnected parts.

`METIS_OPTION_CONTIG` requires the *input* graph to already be connected; requesting
it anyway on a disconnected graph is undefined enough to be unsafe across platforms
-- some METIS builds return an error, others have been observed to crash natively
(a segfault, which no Julia `try`/`catch` can intercept). So contiguity is only ever
requested when `G` is actually connected; a disconnected `G` skips straight to
`relaxedmetisoptions`. When contiguity is requested, weighted partitioning runs
first and retries without vertex weights (`fallbackmetisoptions`) if METIS returns a
single part or throws. Optionally, each resulting part can be split into connected
components.
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

    contiguous = Graphs.is_connected(G)

    part = if contiguous
        try
            partition(
                metisG,
                numberofdivisions;
                options=metisoptions,
                vertexweights=true,
                alg=alg,
            )
        catch e
            e isa MetisError || rethrow()
            nothing
        end
    else
        nothing
    end

    if contiguous && (part === nothing || maximum(part) == minimum(part))
        part === nothing ||
            @warn "METIS partitioning failed to produce multiple parts, retrying without vertex weights"
        part = try
            partition(
                metisG,
                numberofdivisions;
                options=fallbackmetisoptions,
                vertexweights=false,
            )
        catch e
            e isa MetisError || rethrow()
            nothing
        end
    end

    if part === nothing
        contiguous &&
            @warn "METIS could not produce a contiguous partition for a non-contiguous input graph, retrying without the contiguity constraint"
        part = partition(
            metisG,
            numberofdivisions;
            options=relaxedmetisoptions,
            vertexweights=false,
            alg=alg,
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
