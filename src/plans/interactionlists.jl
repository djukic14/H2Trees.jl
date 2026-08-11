# Interaction-list construction for the translating-plan builder.
#
# The generic `WellSeparatedIterator` scans every node at a level. The interaction-list path
# restricts that scan to children of near parents.
#
# The candidate set can be shrunk without changing the answer. `iswellseparated(n, m)` already
# requires `isnear(parent(n), parent(m))`, so every well-separated `m` is a child of a near
# neighbour of `parent(n)`:
#
#     wellsep(n) subseteq union { children(p') : p' in near(parent(n)) }
#
# The predicate is still applied to every surviving candidate, so the result is unchanged.
#
# The near lists this needs are produced by the same top-down sweep, each level from the one
# above, which keeps the whole construction O(N):
#
#     near(root) = {root}
#     near(n)    = { m in children(near(parent(n))) : isnear(n, m) }
#
# This requires near children to imply near parents. Unsupported tree shapes or predicates fall
# back to the original full-level scan.

"""
    _nearpredicate(x)

Return the tree-bound near predicate underlying a translating-nodes iterator, or `nothing`
when it cannot be recovered.

`nothing` is not a failure: it selects the original full-level scan in
[`_buildtranslateplan`](@ref), which is correct for any predicate. Only the standard
well-separated iterators expose a near predicate whose level-monotonicity is established,
so only those take the interaction-list path.
"""
_nearpredicate(::Any) = nothing
function _nearpredicate(functor::_TranslatingFunctor)
    return _nearpredicate(functor.translatingnodesiterator)
end
function _nearpredicate(predicate::_IsTranslatingNodeSingleTreePredicate)
    return _nearpredicate(predicate.TranslatingNodesIterator)
end
function _nearpredicate(builder::_WellSeparatedSingleTreeBuilder)
    return _nearpredicate(builder.iswellseparated)
end
_nearpredicate(predicate::_IsWellSeparatedSingleTreePredicate) = predicate.isnear

"""
    _wellseparatedpredicate(builder)

Return the same-level predicate `builder` applies to its candidates, so the interaction-list
path can call the very same function the full scan would have.
"""
function _wellseparatedpredicate(functor::_TranslatingFunctor)
    return _wellseparatedpredicate(functor.translatingnodesiterator)
end
function _wellseparatedpredicate(predicate::_IsTranslatingNodeSingleTreePredicate)
    return _wellseparatedpredicate(predicate.TranslatingNodesIterator)
end
_wellseparatedpredicate(builder::_WellSeparatedSingleTreeBuilder) = builder.iswellseparated

"""
    _translatinglists(tree, TranslatingNodesIterator)

Return `lists`, a vector indexed by `node - (root(tree) - 1)` holding each node's translating
nodes in ascending node order, or `nothing` when the fast path does not apply.

The per-node lists are exactly what `collect(Int, TranslatingNodesIterator(node))` produces;
only the candidate set the predicate is applied to differs. Returning `nothing` leaves the
caller on its original path.
"""
function _translatinglists(tree, TranslatingNodesIterator)
    supportsnearlists(tree) || return nothing
    isnearpredicate = _nearpredicate(TranslatingNodesIterator)
    isnothing(isnearpredicate) && return nothing
    return _translatinglists(tree, TranslatingNodesIterator, isnearpredicate, Val(true))
end

function _translatinglists(
    tree, TranslatingNodesIterator, isnearpredicate, ::Val{collectfar}, nearout=nothing
) where {collectfar}
    rootoffset = H2Trees.root(tree) - 1
    nnodes = numberofnodes(tree)
    # `nearout` lets a caller keep the near lists this sweep builds anyway (see `_nearlists`),
    # in which case they must not be released level by level.
    nearlists = isnothing(nearout) ? Vector{Vector{Int}}(undef, nnodes) : nearout
    # Only one of these is populated: callers that just need "does this node translate?"
    # (`_translatingflags`) would otherwise pay for materializing every list to then discard it.
    farlists =
        collectfar === true ? Vector{Vector{Int}}(undef, nnodes) : Vector{Vector{Int}}()
    farflags = collectfar === false ? zeros(Bool, nnodes) : Bool[]

    treelevels = collect(H2Trees.levels(tree))
    isempty(treelevels) && return collectfar === true ? farlists : farflags

    # `:nearonly` callers want just the near lists, so no well-separated predicate is consulted.
    iswellseparatedpredicate = if collectfar === :nearonly
        nothing
    else
        _wellseparatedpredicate(TranslatingNodesIterator)
    end

    # Per-node lists are grown by `push!` and then kept, so building them directly would pay
    # doubling overhead on every one: a list ending at ~160 entries passes through
    # 1+2+4+...+256 = 511 slots. Measured on 400k points that was 93.6MiB allocated to store a
    # 22.9MiB payload. Filling a reused scratch buffer and copying it out once costs exactly the
    # payload instead.
    candidatebuffer = Int[]
    nearscratch = Int[]
    farscratch = Int[]
    for (levelid, treelevel) in enumerate(treelevels)
        levelnodes = LevelIterator(tree, treelevel)
        isempty(levelnodes) && continue

        # A node's candidates come from its parent's near list, so once a level is done every
        # near list two levels up is dead. Releasing them keeps peak memory at ~2 levels rather
        # than the whole tree.
        if levelid > 2 && isnothing(nearout)
            for stale in LevelIterator(tree, treelevels[levelid - 2])
                nearlists[stale - rootoffset] = EMPTYNODELIST
            end
        end

        for node in levelnodes
            # The coarsest level has no parent level to restrict against. It is also tiny (one
            # node for a full tree), so scanning it outright costs nothing.
            candidates = if levelid == 1
                levelnodes
            else
                _appendchildrenofnearparents!(
                    candidatebuffer, tree, node, nearlists, rootoffset
                )
            end

            empty!(nearscratch)
            empty!(farscratch)
            for candidate in candidates
                # Check both argument orders so user predicates need not be symmetric. If the
                # sweep predicate ever differs from the well-separated predicate's `isnear`,
                # split this short-circuit into separate candidate and rejection tests.
                if isnearpredicate(tree, node, candidate) ||
                    isnearpredicate(tree, candidate, node)
                    push!(nearscratch, candidate)
                elseif collectfar !== :nearonly &&
                    iswellseparatedpredicate(tree, candidate, node)
                    if collectfar === true
                        push!(farscratch, candidate)
                    else
                        farflags[node - rootoffset] = true
                    end
                end
            end

            # `LevelIterator` yields ascending node ids, so the original lists are ascending.
            # The candidate set here is grouped by parent instead, so restore that order.
            # `copy` sizes each stored list exactly, unlike growing it by `push!`.
            nearlists[node - rootoffset] = copy(sort!(nearscratch))
            collectfar === true && (farlists[node - rootoffset] = copy(sort!(farscratch)))
        end
    end

    return collectfar === true ? farlists : farflags
end

# Shared placeholder for released near lists.
const EMPTYNODELIST = Int[]

"""
    _nearlists(tree, isnear)

Return per-node near-node lists indexed by `node - (root(tree) - 1)`, in ascending node order,
or `nothing` when no near predicate can be recovered from `isnear`.

Same top-down sweep as [`_translatinglists`](@ref), but retaining near lists for every level.
"""
function _nearlists(tree, isnearpredicate)
    nearout = Vector{Vector{Int}}(undef, numberofnodes(tree))
    # `Val(:nearonly)` skips the well-separated branch entirely, so no such predicate is needed.
    _translatinglists(tree, nothing, isnearpredicate, Val(:nearonly), nearout)
    return nearout
end

"""
    _translatingflagsfromlists(lists)

Derive the per-node "does this node translate?" flags from already-built interaction lists,
so a caller building both plans of a pair does not sweep the tree twice. Propagates `nothing`.
"""
_translatingflagsfromlists(::Nothing) = nothing
_translatingflagsfromlists(lists) = [!isempty(list) for list in lists]

"""
    _sharesiterator(aggregatenode, translatingnodesiterator)

Return whether `aggregatenode` is `istranslatingnode` over exactly `translatingnodesiterator`,
so that "does this node translate?" and the interaction lists answer the same question.

This must be checked, not assumed. `aggregatenode` is a free choice: `aggregateallnodes` stores
every node and `aggregaterootonly` stores just the root, and substituting translating-flags for
either silently builds the wrong aggregation plan. Only the default `PlanBuilder` wiring, which
threads one iterator into both, may share.
"""
_sharesiterator(::Any, ::Any) = false
function _sharesiterator(
    aggregatenode::_IsTranslatingNodePredicateBuilder, translatingnodesiterator
)
    return aggregatenode.TranslatingNodesIterator === translatingnodesiterator
end

"""
    _sharedtranslatingflags(aggregatenode, translatingnodesiterator, lists)

Flags derived from `lists` when [`_sharesiterator`](@ref) permits it, otherwise `nothing` so the
plan builder evaluates the caller's own `aggregatenode` predicate.
"""
function _sharedtranslatingflags(aggregatenode, translatingnodesiterator, lists)
    _sharesiterator(aggregatenode, translatingnodesiterator) || return nothing
    return _translatingflagsfromlists(lists)
end

"""
    _translatingflags(tree, storenode)

Return a `Bool` vector indexed by `node - (root(tree) - 1)`, true where `storenode(node)` is,
or `nothing` when the fast path does not apply.

`istranslatingnode` answers "does this node translate to anything?" by scanning the level until
the first hit. Deriving the answer from interaction lists makes it a lookup. Prefer
[`_translatingflagsfromlists`](@ref) when the lists are already at hand.
"""
function _translatingflags(tree, storenode)
    supportsnearlists(tree) || return nothing
    isnearpredicate = _nearpredicate(storenode)
    isnothing(isnearpredicate) && return nothing
    return _translatinglists(tree, storenode, isnearpredicate, Val(false))
end

"""
    _appendchildrenofnearparents!(buffer, tree, node, nearlists, rootoffset)

Fill `buffer` with the children of every near neighbour of `node`'s parent and return it.

This is the candidate superset of `node`'s well-separated nodes; `buffer` is reused across
nodes, so the result is only valid until the next call.
"""
function _appendchildrenofnearparents!(buffer, tree, node, nearlists, rootoffset)
    empty!(buffer)
    nodeparent = H2Trees.parent(tree, node)
    iszero(nodeparent) && return buffer

    for nearparent in nearlists[nodeparent - rootoffset]
        for child in H2Trees.children(tree, nearparent)
            push!(buffer, child)
        end
    end
    return buffer
end
