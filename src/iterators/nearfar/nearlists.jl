
# Cached near lists.
#
# `NearNodeIterator` normally scans a whole level. `_nearlists` builds candidate sets in one
# top-down sweep:
#
#     near(root) = {root}
#     near(n)    = { m in children(near(parent(n))) : isnear(n, m) or isnear(m, n) }
#
# The cached list is only a candidate set; the caller's predicate is still applied.
# Two properties are required:
#
#   1. The sweep predicate must be level-monotone: a near child implies near parents.
#      Unknown predicates return `nothing`, selecting the original scan.
#   2. A consumer's predicate must be no wider than the one the cache was built from.
#      `_LeafPredicateFunctor` narrows (`isleaf && isnear`), so it may filter a cache built
#      from its inner predicate, but must not seed a sweep itself.
#
# Monotonicity is a property of the TREE SHAPE as much as of the predicate, which is why
# `supportsnearlists` also gates on the tree trait.

"""
    nearcandidatepredicate(predicate)

Return a level-monotone near predicate whose near set contains `predicate`'s, or `nothing`
when no such predicate is known for it.

`nothing` is not a failure: it selects the original full-level scan, which is correct for any
predicate. A non-`nothing` result means the returned predicate contains the original one. A
wrapper that only *narrows* its inner predicate (as `_LeafPredicateFunctor` does) satisfies
this by unwrapping to that inner predicate, which is the extension point for predicates
defined outside H2Trees.

Whether the returned predicate is also level-monotone, the other half of what the sweep
needs, depends on the tree's shape rather than on the predicate alone and is decided
separately by [`supportsnearlists`](@ref). The library's geometric predicates are recognised
here for that reason: on the shapes `supportsnearlists` admits they are monotone, and on the
shapes it rejects they never get the chance.
"""
nearcandidatepredicate(::Any) = nothing
nearcandidatepredicate(predicate::typeof(isnear)) = predicate
nearcandidatepredicate(predicate::IsNearNotBlockTreeFunctor) = predicate
function nearcandidatepredicate(predicate::_LeafPredicateFunctor)
    return nearcandidatepredicate(predicate.predicate)
end

"""
    supportsnearlists(tree)

Return whether the top-down near recursion is sound for `tree`'s shape.

Only box trees ([`isTwoNTree`](@ref)) qualify. There, a parent's box contains its children's
and the axis-aligned gap can only shrink going up while [`isneargap`](@ref)'s allowance grows
with the half-size, so a near child forces near parents.

Ball trees keep the original scan. Their current radius predicate is monotone under
containment, but the measured speedup was small enough that enabling the cache remains a
separate performance and accuracy decision.

A shape not covered here keeps the original full-level scan, which is correct for any tree.
"""
supportsnearlists(tree) = supportsnearlists(treetrait(tree))
supportsnearlists(::AbstractTreeTrait) = false
supportsnearlists(::isTwoNTree) = true

"""
    NearListCache

Per-node near-node lists, indexed by node id, used as the candidate set for near iteration.

Built by [`nearlistcache`](@ref). Each entry is a superset of the node's true near nodes, so
the consuming predicate must still be applied; see the discussion in `nearlists.jl`.
"""
struct NearListCache
    lists::Vector{Vector{Int}}
    rootoffset::Int
end

Base.getindex(cache::NearListCache, node::Int) = cache.lists[node - cache.rootoffset]

"""
    nearlistcache(tree, predicate)

Return a [`NearListCache`](@ref) usable as the candidate set for `predicate`, or `nothing`
when the full-level scan must be kept.

`nothing` is returned for a tree shape [`supportsnearlists`](@ref) rejects (which includes
block trees, whose pairs span two independently built trees a single-tree sweep cannot
enumerate), for a predicate [`nearcandidatepredicate`](@ref) does not vouch for, and for a
tree whose nodes are not all covered by `levels(tree)`.
"""
function nearlistcache(tree, predicate)
    supportsnearlists(tree) || return nothing

    candidatepredicate = nearcandidatepredicate(predicate)
    isnothing(candidatepredicate) && return nothing

    lists = _nearlists(tree, candidatepredicate)
    # The sweep only visits nodes reachable through `levels(tree)`. A node outside that has no
    # list, and an unfilled slot must send the caller back to the scan rather than read as an
    # empty near set, which would silently drop interactions.
    all(i -> isassigned(lists, i), eachindex(lists)) || return nothing

    return NearListCache(lists, root(tree) - 1)
end

"""
    _nearcandidates(nearlists, tree, node)

The set of nodes a near filter for `node` should be applied to: every node at its level when
`nearlists` is `nothing`, or the cached candidate list otherwise.
"""
_nearcandidates(::Nothing, tree, node) = LevelIterator(tree, level(tree, node))
_nearcandidates(cache::NearListCache, tree, node) = cache[node]
