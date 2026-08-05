# Shared construction algorithm for `AggregatePlan` and `DisaggregatePlan`.
#
# Both plans visit nodes selected by a storage predicate plus descendants whose
# ancestors are selected. The differences are traversal order, level indexing,
# level validation, empty-plan policy, and the final concrete plan type.
abstract type StoredNodePlanBuildKind end
struct AggregateStoredNodeBuild <: StoredNodePlanBuildKind end
struct DisaggregateStoredNodeBuild <: StoredNodePlanBuildKind end

_storednodebuildlevels(::AggregateStoredNodeBuild, tree) = collect(H2Trees.levels(tree))
function _storednodebuildlevels(::DisaggregateStoredNodeBuild, tree)
    return collect(reverse(H2Trees.levels(tree)))
end

function _storednodelevelid(::AggregateStoredNodeBuild, tree, level)
    return numberoflevels(tree) - level + minimumlevel(tree)
end
_storednodelevelid(::DisaggregateStoredNodeBuild, tree, level) = leveltolevelid(tree, level)

function _validatestorednodelevels(::AggregateStoredNodeBuild, levels)
    return _validatedaggregationlevels(levels)
end
function _validatestorednodelevels(::DisaggregateStoredNodeBuild, levels)
    return _validateddisaggregationlevels(levels)
end

function _checkemptystorednodeplan(::AggregateStoredNodeBuild, nodes, levels)
    (isempty(nodes) || isempty(levels)) &&
        throw(ArgumentError("Empty AggregatePlan not supported."))
    return nothing
end
_checkemptystorednodeplan(::DisaggregateStoredNodeBuild, nodes, levels) = nothing

function _materializestorednodeplan(
    ::AggregateStoredNodeBuild, nodes, levels, storenodes, rootoffset, tree
)
    return AggregatePlan(nodes, levels, storenodes, rootoffset, tree)
end
function _materializestorednodeplan(
    ::DisaggregateStoredNodeBuild, nodes, levels, storenodes, rootoffset, tree
)
    return DisaggregatePlan(nodes, levels, storenodes, rootoffset, tree)
end

function _buildstorednodeplan(kind::StoredNodePlanBuildKind, tree, storenode)
    storenodesarray = zeros(Bool, numberofnodes(tree))
    root = H2Trees.root(tree)
    rootoffset = root - 1

    rawlevels = zeros(Int, numberoflevels(tree))
    visitnodes = Vector{Vector{Int}}(undef, numberoflevels(tree))

    lk = Threads.SpinLock()
    for level in _storednodebuildlevels(kind, tree)
        levelnodes = Int[]

        @threads for node in LevelIterator(tree, level)
            nodeindex = node - rootoffset

            if storenode(node)
                storenodesarray[nodeindex] = true
                @lock lk push!(levelnodes, node)
            end

            storenodesarray[nodeindex] && continue

            for parent in ParentUpwardsIterator(tree, node)
                if storenode(parent)
                    @lock lk push!(levelnodes, node)
                    break
                end
            end
        end

        isempty(levelnodes) && continue

        levelid = _storednodelevelid(kind, tree, level)
        rawlevels[levelid] = level
        visitnodes[levelid] = levelnodes
    end

    indicestodelete = Int[]
    for i in eachindex(visitnodes)
        isassigned(visitnodes, i) || push!(indicestodelete, i)
    end
    deleteat!(rawlevels, indicestodelete)
    deleteat!(visitnodes, indicestodelete)

    rawlevels = _validatestorednodelevels(kind, rawlevels)
    _checkemptystorednodeplan(kind, visitnodes, rawlevels)

    return _materializestorednodeplan(
        kind, visitnodes, rawlevels, storenodesarray, rootoffset, tree
    )
end
