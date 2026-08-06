"""
    TranslatePlanBuildKind

Internal strategy marker for the shared translating-plan builder.

`_buildtranslateplan` owns the traversal used by both
[`AggregateTranslatePlan`](@ref) and [`DisaggregateTranslatePlan`](@ref). Each
build kind supplies only the plan-specific hooks: level storage order, dictionary
recording, node-flag construction, empty-plan policy, and final materialization.
"""
abstract type TranslatePlanBuildKind end

"""
    AggregateTranslateBuild <: TranslatePlanBuildKind

Build-kind marker for [`AggregateTranslatePlan`](@ref).
"""
struct AggregateTranslateBuild <: TranslatePlanBuildKind end

"""
    DisaggregateTranslateBuild <: TranslatePlanBuildKind

Build-kind marker for [`DisaggregateTranslatePlan`](@ref).
"""
struct DisaggregateTranslateBuild <: TranslatePlanBuildKind end

_translatebuildlevels(::TranslatePlanBuildKind, tree) = collect(H2Trees.levels(tree))

function _translatelevelid(::AggregateTranslateBuild, tree, level)
    return numberoflevels(tree) - level + minimumlevel(tree)
end
_translatelevelid(::DisaggregateTranslateBuild, tree, level) = leveltolevelid(tree, level)

function _validatetranslatelevels(::AggregateTranslateBuild, levels)
    return _validatedaggregationlevels(levels)
end
function _validatetranslatelevels(::DisaggregateTranslateBuild, levels)
    return _validateddisaggregationlevels(levels)
end

"""
    _recordtranslatenode!(kind, leveldict, nodeflags, node, nodeindex, tfnodes)

Record one visited translating-plan node into the per-level dictionary and node
flags for `kind`.
"""
function _recordtranslatenode!(
    ::AggregateTranslateBuild, leveldict, nodeflags, node, nodeindex, tfnodes
)
    for tfnode in tfnodes
        if !haskey(leveldict, tfnode)
            leveldict[tfnode] = [node]
        else
            push!(leveldict[tfnode], node)
        end
    end
    return nothing
end

function _recordtranslatenode!(
    ::DisaggregateTranslateBuild, leveldict, nodeflags, node, nodeindex, tfnodes
)
    isempty(tfnodes) || (leveldict[nodeindex] = tfnodes)
    nodeflags[nodeindex] = true
    return nothing
end

_initialnodeflags(::AggregateTranslateBuild, tree) = nothing
_initialnodeflags(::DisaggregateTranslateBuild, tree) = zeros(Bool, numberofnodes(tree))

"""
    _finalnodeflags(kind, nodeflags, dicts, tree)

Return the node-flag vector stored by the materialized translating plan.
"""
function _finalnodeflags(::AggregateTranslateBuild, nodeflags, dicts, tree)
    return _computeistranslatingnodes(dicts, tree)
end
_finalnodeflags(::DisaggregateTranslateBuild, nodeflags, dicts, tree) = nodeflags

function _checkemptytranslateplan(::AggregateTranslateBuild, nodes, levels)
    (isempty(nodes) || isempty(levels)) &&
        throw(ArgumentError("Empty AggregatePlan not supported."))
    return nothing
end
_checkemptytranslateplan(::DisaggregateTranslateBuild, nodes, levels) = nothing

function _materializetranslateplan(
    ::AggregateTranslateBuild, dicts, nodes, levels, flags, rootoffset, tree
)
    return _makeaggregatetranslateplan(dicts, nodes, levels, flags, rootoffset, tree)
end
function _materializetranslateplan(
    ::DisaggregateTranslateBuild, dicts, nodes, levels, flags, rootoffset, tree
)
    return _makedisaggregatetranslateplan(dicts, nodes, levels, flags, rootoffset, tree)
end

"""
    _buildtranslateplan(kind::TranslatePlanBuildKind, tree, TranslatingNodesIterator)

Shared construction algorithm for `AggregateTranslatePlan` (`kind = AggregateTranslateBuild()`)
and `DisaggregateTranslatePlan` (`kind = DisaggregateTranslateBuild()`). Visits every tree level,
marks a node for visiting if it (or an ancestor) has a nonempty `TranslatingNodesIterator`, and
records the per-level translating/receiving pairs into a compacted, validated set of levels.

The builder preserves each concrete plan's storage convention: aggregation
levels are stored from fine to coarse, while disaggregation levels are stored
from coarse to fine.
"""
function _buildtranslateplan(kind::TranslatePlanBuildKind, tree, TranslatingNodesIterator)
    levels = _translatebuildlevels(kind, tree)
    rootoffset = H2Trees.root(tree) - 1

    rawlevels = zeros(Int, numberoflevels(tree))
    visitnodes = Vector{Vector{Int}}(undef, numberoflevels(tree))
    dicts = Vector{Dict{Int,Vector{Int}}}(undef, numberoflevels(tree))
    nodeflags = _initialnodeflags(kind, tree)

    lk = Threads.SpinLock()
    for level in levels
        levelnodes = Int[]
        leveldict = Dict{Int,Vector{Int}}()

        @threads for node in LevelIterator(tree, level)
            nodeindex = node - rootoffset
            nodehastobevisited = false

            tfnodes = collect(Int, TranslatingNodesIterator(node))
            !isempty(tfnodes) && (nodehastobevisited = true)

            if !nodehastobevisited
                for parent in ParentUpwardsIterator(tree, node)
                    for _ in TranslatingNodesIterator(parent)
                        nodehastobevisited = true
                        break
                    end
                end
            end

            if nodehastobevisited
                lock(lk) do
                    push!(levelnodes, node)
                    return _recordtranslatenode!(
                        kind, leveldict, nodeflags, node, nodeindex, tfnodes
                    )
                end
            end
        end

        levelid = _translatelevelid(kind, tree, level)
        (isempty(levelnodes) && isempty(leveldict)) && continue
        rawlevels[levelid] = level
        visitnodes[levelid] = levelnodes
        dicts[levelid] = leveldict
    end

    indicestodelete = Int[]
    for i in eachindex(visitnodes)
        isassigned(visitnodes, i) || push!(indicestodelete, i)
    end
    deleteat!(rawlevels, indicestodelete)
    deleteat!(visitnodes, indicestodelete)
    deleteat!(dicts, indicestodelete)

    rawlevels = _validatetranslatelevels(kind, rawlevels)
    _checkemptytranslateplan(kind, visitnodes, rawlevels)

    nodeflags = _finalnodeflags(kind, nodeflags, dicts, tree)

    return _materializetranslateplan(
        kind, dicts, visitnodes, rawlevels, nodeflags, rootoffset, tree
    )
end
