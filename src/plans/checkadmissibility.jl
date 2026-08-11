
# Admissibility checking ###################################################################
#
# Validate assembled translation plans against near/far geometry, not against a matvec.

"""
    AdmissibilityFinding

One problem reported by [`checkadmissibility`](@ref).

`severity` is `:error` for a plan that contradicts its own `isnear` predicate or references
nodes that cannot be resolved, and `:warning` for a plan that is internally consistent but
contains geometrically marginal far pairs.

`kind` is one of `:nearpairtranslated`, `:marginalgap`, `:wronglevel`, `:unknownnode`,
`:coveragegap`, `:coverageduplicate`.

Gap fields are `missing` for findings that are not about geometry.
"""
struct AdmissibilityFinding
    severity::Symbol
    kind::Symbol
    level::Union{Int,Missing}
    receivingnode::Union{Int,Missing}
    translatingnode::Union{Int,Missing}
    gap::Union{Float64,Missing}
    gapoverminnodesize::Union{Float64,Missing}
    gapovermaxnodesize::Union{Float64,Missing}
    message::String
end

"""
    AdmissibilityReport

Result of [`checkadmissibility`](@ref). `ok` is true when no `:error` finding was produced.
`:warning` findings do not clear it.

Inspect `findings` directly, or let `show` summarise it.
"""
struct AdmissibilityReport
    ok::Bool
    findings::Vector{AdmissibilityFinding}
end

"""
    haserrors(report)

Return whether an [`AdmissibilityReport`](@ref) contains at least one `:error`
finding.
"""
haserrors(r::AdmissibilityReport) = any(f -> f.severity === :error, r.findings)

"""
    haswarnings(report)

Return whether an [`AdmissibilityReport`](@ref) contains at least one `:warning`
finding.
"""
haswarnings(r::AdmissibilityReport) = any(f -> f.severity === :warning, r.findings)

function Base.show(io::IO, ::MIME"text/plain", r::AdmissibilityReport)
    nerr = count(f -> f.severity === :error, r.findings)
    nwarn = count(f -> f.severity === :warning, r.findings)
    println(io, "AdmissibilityReport  ok=", r.ok, "  errors=", nerr, "  warnings=", nwarn)
    for f in r.findings
        print(io, "  [", f.severity, "] ", f.kind, ": ", f.message)
        if !ismissing(f.gap)
            print(io, "  (gap=", round(f.gap; sigdigits=3))
            ismissing(f.gapoverminnodesize) ||
                print(io, ", gap/min(size)=", round(f.gapoverminnodesize; sigdigits=3))
            ismissing(f.gapovermaxnodesize) ||
                print(io, ", gap/max(size)=", round(f.gapovermaxnodesize; sigdigits=3))
            print(io, ")")
        end
        println(io)
    end
    return nothing
end

"""
    boxgap(center_a, center_b, halfsize_a, halfsize_b)

Axis-aligned gap between two boxes: the Euclidean length of the per-axis overhangs, counting
only axes where the boxes do not overlap. Zero when they touch or intersect.

This is the raw geometric quantity the admissibility predicate is built on (see
[`isneargap`](@ref)) exposed separately so diagnostics can report the real gap independently
of whatever predicate a plan happened to be built with.
"""
function boxgap(
    center_a::AbstractVector, center_b::AbstractVector, halfsize_a::T, halfsize_b
) where {T}
    reach = halfsize_a + halfsize_b
    gapsquared = zero(T)
    for i in eachindex(center_a)
        gap = abs(center_a[i] - center_b[i]) - reach
        gap > 0 && (gapsquared += gap * gap)
    end
    return sqrt(gapsquared)
end

"""
    ballgap(center_a, center_b, radius_a, radius_b)

Separation between two bounding balls: centre-to-centre distance minus both radii, floored at
zero when the balls touch or overlap.

The [`BoundingBallTree`](@ref)/[`KMeansTree`](@ref) counterpart to [`boxgap`](@ref): see
[`nodegap`](@ref) for the dispatch that picks between the two.
"""
function ballgap(
    center_a::AbstractVector, center_b::AbstractVector, radius_a::T, radius_b
) where {T}
    return max(zero(T), norm(center_a - center_b) - radius_a - radius_b)
end

"""
    nodesize(tree, node)

The characteristic bounding-volume size `checkadmissibility`'s gap margin is measured in:
`halfsize` for a box-shaped tree ([`isTwoNTree`](@ref)), `radius` for a ball-shaped one
([`isBoundingBallTree`](@ref)).
"""
nodesize(tree, node) = nodesize(tree, node, treetrait(tree))
nodesize(tree, node, ::isTwoNTree) = halfsize(tree, node)
nodesize(tree, node, ::isBoundingBallTree) = radius(tree, node)

"""
    nodegap(treea, nodea, treeb, nodeb)

Geometric separation between two tree nodes, dispatched on `treea`'s own shape: [`boxgap`](@ref)
for a box-shaped tree, [`ballgap`](@ref) for a ball-shaped one. `treea` and `treeb` are always the
same shape in practice (a `BlockTree`'s two sides are both `TwoNTree`s by construction) so
dispatching on `treea` alone is sufficient.
"""
function nodegap(treea, nodea, treeb, nodeb)
    return nodegap(treea, nodea, treeb, nodeb, treetrait(treea))
end
function nodegap(treea, nodea, treeb, nodeb, ::isTwoNTree)
    return boxgap(
        center(treea, nodea),
        center(treeb, nodeb),
        halfsize(treea, nodea),
        halfsize(treeb, nodeb),
    )
end
function nodegap(treea, nodea, treeb, nodeb, ::isBoundingBallTree)
    return ballgap(
        center(treea, nodea),
        center(treeb, nodeb),
        radius(treea, nodea),
        radius(treeb, nodeb),
    )
end

# `isnear` may be passed in two shapes:
#
#   - UNRESOLVED: a factory taking the tree and returning the predicate. `H2Trees.isnear()`
#     returns an `IsNearFunctor` of this kind, and a user closure `tree -> predicate` is the same
#     shape.
#   - ALREADY RESOLVED: the node-comparison callable itself. `H2Trees.isnear(A)` on an H2 map
#     hands back the concrete `IsNearBlockTreeFunctor`/`IsNearNotBlockTreeFunctor` the map was built
#     with, and a hand-written `(tree, a, b) -> ...` / `(testtree, trialtree, a, b) -> ...` closure
#     is equally natural to pass to a diagnostic.
#
# Telling them apart by "can it be called on the tree alone?" covers both the library types and
# arbitrary user closures, which typing the two known functors could not. The named methods below
# are kept anyway: they are the common path, and they document the intent at a glance.
function _nearpredicate(isnearfunctor, tree)
    applicable(isnearfunctor, tree) || return isnearfunctor
    return isnearfunctor(tree)
end

_nearpredicate(f::IsNearBlockTreeFunctor, tree) = f
_nearpredicate(f::IsNearNotBlockTreeFunctor, tree) = f

# Validated ONCE, up front, rather than per pair: a shape mismatch otherwise surfaces as a
# `MethodError` thrown from deep inside the translation-pair loop, where the reported argument types
# say nothing about which of the accepted forms the caller should have used.
function _assertpredicateshape(pred, tree, trait)
    _predicateapplicable(pred, tree, trait) && return nothing
    return error(
        "the `isnear` predicate cannot be called with this tree's node-comparison signature " *
        "($(trait isa isBlockTree ? "isnear(testtree, trialtree, testnode, trialnode)" :
            "isnear(tree, nodea, nodeb)")). Accepted forms: a factory taking the tree and " *
        "returning the predicate (e.g. `H2Trees.isnear()`, or `tree -> predicate`), or the " *
        "resolved predicate itself (e.g. `H2Trees.isnear(A)` from a map, or a plain closure of " *
        "the signature above). Got a $(typeof(pred)).",
    )
end

function _predicateapplicable(pred, tree, ::AbstractTreeTrait)
    return applicable(pred, tree, 1, 1)
end

function _predicateapplicable(pred, tree, ::isBlockTree)
    return applicable(pred, testtree(tree), trialtree(tree), 1, 1)
end

function _callnear(
    pred, tree, receivingtree, translatingtree, rnode, tnode, ::AbstractTreeTrait
)
    return pred(receivingtree, rnode, tnode)
end

function _callnear(pred, tree, receivingtree, translatingtree, rnode, tnode, ::isBlockTree)
    if receivingtree === testtree(tree)
        return pred(receivingtree, translatingtree, rnode, tnode)
    else
        return pred(translatingtree, receivingtree, tnode, rnode)
    end
end

"""
    checkadmissibility(tree, plans; isnear=H2Trees.isnear(), mingapboxes=DEFAULTNEARGAPBOXES,
                       coverage=true, throw=true, maxfindings=200)

Validate an assembled translation scheme against its near/far geometry.

`plans` is normally the [`PlanSet`](@ref) returned by [`buildplans`](@ref).
`isnear` may be either a factory such as `H2Trees.isnear()` or an already
resolved node-comparison predicate.

The report contains `:error` findings for impossible or inconsistent plan data
(near pairs scheduled as far, wrong levels, unknown nodes, coverage gaps or
duplicates) and `:warning` findings for marginal far pairs whose geometric gap
is below `mingapboxes * min(nodesize)`.

With `coverage=true`, the near field is compared with the far values actually
scheduled by the translating plan, including translations addressed to ancestors
and carried down by disaggregation. This is intentionally stronger than
recomputing the far field from `FarNodeIterator`, but it is also expensive; use
`coverage=false` on large diagnostic runs.

`throw=true` raises on errors and warns on warnings. Use `throw=false` to inspect
the returned [`AdmissibilityReport`](@ref). This is a diagnostic, not a matvec
hot-path check.
"""
function checkadmissibility(
    tree,
    plans;
    isnear=isnear(),
    mingapboxes=DEFAULTNEARGAPBOXES,
    coverage::Bool=true,
    throw::Bool=true,
    maxfindings::Int=200,
)
    findings = AdmissibilityFinding[]
    findinglock = ReentrantLock()
    trait = treetrait(tree)
    pred = _nearpredicate(isnear, tree)
    _assertpredicateshape(pred, tree, trait)

    for plan in _translatingplans(plans)
        _checkplanadmissibility!(
            findings, findinglock, tree, plan, pred, trait, mingapboxes, maxfindings
        )
        if coverage && length(findings) < maxfindings
            _checkplancoverage!(findings, findinglock, tree, plan, pred, trait, maxfindings)
        end
        length(findings) >= maxfindings && break
    end

    report = AdmissibilityReport(!any(f -> f.severity === :error, findings), findings)
    throw && _raise(report)
    return report
end

# Both traversal directions carry a translating plan (forward and its adjoint counterpart);
# check whichever ones are present. `translatingplan` already errors if a pair does not have
# exactly one translating member, so this also validates the pairing itself.
"""
    _translatingplans(plans)

Return the translating plans exposed by a plan set.

Both the forward pair and the adjoint pair are checked when present. The helper
also validates that each pair contains exactly one translating plan.
"""
function _translatingplans(plans)
    out = Any[]
    if hasproperty(plans, :trialaggregationplan) &&
        hasproperty(plans, :testdisaggregationplan)
        push!(
            out, translatingplan(plans.trialaggregationplan, plans.testdisaggregationplan)
        )
    end
    if hasproperty(plans, :testaggregationplan) &&
        hasproperty(plans, :trialdisaggregationplan)
        push!(
            out, translatingplan(plans.testaggregationplan, plans.trialdisaggregationplan)
        )
    end
    isempty(out) && return error(
        "`plans` exposes no translating plan pair: expected `trialaggregationplan` + " *
        "`testdisaggregationplan` (and optionally `testaggregationplan` + " *
        "`trialdisaggregationplan`), as returned by `buildplans`.",
    )
    return out
end

"""
    _checkplanadmissibility!(findings, findinglock, tree, plan, pred, trait,
        mingapboxes, maxfindings)

Check every translation pair scheduled by one translating plan.

Findings are appended for pairs filed at the wrong level, unresolved nodes, near
pairs scheduled as far interactions, and far pairs whose geometric gap is below
the configured margin.
"""
function _checkplanadmissibility!(
    findings, findinglock, tree, plan, pred, trait, mingapboxes, maxfindings
)
    receivingtree_ = receivingtree(tree, plan)
    translatingtree_ = translatingtree(tree, plan)
    jobs = Tuple{Int,Int}[]
    for level in levels(plan)
        for rnode in receivingnodes(plan, level)
            push!(jobs, (level, rnode))
        end
    end

    Threads.@threads for i in eachindex(jobs)
        level, rnode = jobs[i]
        _checkreceivingnodeadmissibility!(
            findings,
            findinglock,
            tree,
            plan,
            pred,
            trait,
            mingapboxes,
            maxfindings,
            receivingtree_,
            translatingtree_,
            level,
            rnode,
        )
    end
    return nothing
end

function _checkreceivingnodeadmissibility!(
    findings,
    findinglock,
    tree,
    plan,
    pred,
    trait,
    mingapboxes,
    maxfindings,
    receivingtree,
    translatingtree,
    level,
    rnode,
)
    _checknodelevel!(
        findings, findinglock, receivingtree, rnode, level, :receiving, maxfindings
    ) || return nothing
    rh = nodesize(receivingtree, rnode)

    for tnode in translatingnodes(plan, rnode, level)
        _findingsfull(findings, findinglock, maxfindings) && return nothing
        _checknodelevel!(
            findings, findinglock, translatingtree, tnode, level, :translating, maxfindings
        ) || continue

        th = nodesize(translatingtree, tnode)
        gap = nodegap(receivingtree, rnode, translatingtree, tnode)
        minh = min(rh, th)
        maxh = max(rh, th)

        if _callnear(pred, tree, receivingtree, translatingtree, rnode, tnode, trait)
            _pushfinding!(
                findings,
                findinglock,
                AdmissibilityFinding(
                    :error,
                    :nearpairtranslated,
                    level,
                    rnode,
                    tnode,
                    Float64(gap),
                    Float64(gap / minh),
                    Float64(gap / maxh),
                    "pair is scheduled for translation but `isnear` reports it as near",
                ),
                maxfindings,
            )
        elseif gap <= mingapboxes * minh
            _pushfinding!(
                findings,
                findinglock,
                AdmissibilityFinding(
                    :warning,
                    :marginalgap,
                    level,
                    rnode,
                    tnode,
                    Float64(gap),
                    Float64(gap / minh),
                    Float64(gap / maxh),
                    "far pair is separated by less than the recommended $(mingapboxes) x min(nodesize)",
                ),
                maxfindings,
            )
        end
    end
    return nothing
end

# Returns false when the node could not be validated, so the caller skips geometry on it rather
# than throwing out of `center`/`halfsize`.
"""
    _checknodelevel!(findings, findinglock, tree, node, planlevel, role, maxfindings)

Validate that `node` exists in `tree` and is filed under its actual level.

Returns `false` after recording a finding when geometric checks for that node
should be skipped.
"""
function _checknodelevel!(findings, findinglock, tree, node, planlevel, role, maxfindings)
    _findingsfull(findings, findinglock, maxfindings) && return false
    actual = try
        level(tree, node)
    catch
        _pushfinding!(
            findings,
            findinglock,
            AdmissibilityFinding(
                :error,
                :unknownnode,
                planlevel,
                role === :receiving ? node : missing,
                role === :translating ? node : missing,
                missing,
                missing,
                missing,
                "$(role) node $(node) cannot be resolved in its tree",
            ),
            maxfindings,
        )
        return false
    end
    if actual != planlevel
        _pushfinding!(
            findings,
            findinglock,
            AdmissibilityFinding(
                :error,
                :wronglevel,
                planlevel,
                role === :receiving ? node : missing,
                role === :translating ? node : missing,
                missing,
                missing,
                missing,
                "$(role) node $(node) is filed at plan level $(planlevel) but sits at tree level $(actual)",
            ),
            maxfindings,
        )
        return false
    end
    return true
end

"""
    _staticchunks(n, nchunks)

Split `1:n` into contiguous, roughly even ranges.

Used by the coverage pass to give each spawned task its own scratch buffer
without relying on `Threads.@threads :static`, which cannot be nested or run
concurrently.
"""
function _staticchunks(n::Int, nchunks::Int)
    nchunks = clamp(nchunks, 1, max(1, n))
    base, extra = divrem(n, nchunks)
    chunks = Vector{UnitRange{Int}}(undef, nchunks)
    lo = 1
    for c in 1:nchunks
        len = base + (c <= extra ? 1 : 0)
        chunks[c] = lo:(lo + len - 1)
        lo += len
    end
    return chunks
end

"""
    _checkplancoverage!(findings, findinglock, tree, plan, nearpred, trait, maxfindings)

Check that near-field values plus plan-scheduled far-field values cover the
translating side exactly once for every receiving leaf.

The far set is read from the plan, including translations scheduled on leaf
ancestors, so omissions or duplicates in the assembled plan are detected
directly.
"""
function _checkplancoverage!(
    findings, findinglock, tree, plan, nearpred, trait, maxfindings
)
    receivingtree_ = receivingtree(tree, plan)
    translatingtree_ = translatingtree(tree, plan)
    planlevels = levels(plan)
    expected = _allvalues(translatingtree_)
    leafnodes = collect(leaves(receivingtree_))
    # Built once per plan rather than rescanning the receiving level for every leaf; `nothing`
    # (block trees, unrecognised predicates) keeps the original scan. See `nearlistcache`.
    nearlists = nearlistcache(receivingtree_, nearpred)
    # `nothing` when the value ids are too sparse or not unique; then each leaf sorts instead.
    coverageindex = _coverageindex(expected)

    # One `got` buffer per chunk, reused across every leaf in that chunk, instead of a fresh
    # allocation per leaf; each chunk's task owns its buffer for the task's whole lifetime, so
    # no lock is needed for `got` itself (the findings lock is separate and unchanged).
    @sync for chunk in _staticchunks(length(leafnodes), Threads.nthreads())
        Threads.@spawn begin
            got = Int[]
            sizehint!(got, length(expected))
            # Same reasoning as `got`: the stamp array is this task's private scratch, reused
            # across its leaves. `generation` makes reuse safe without clearing it:
            # a slot counts as covered only when it carries the current leaf's number.
            stamps = isnothing(coverageindex) ? Int[] : zeros(Int, _nslots(coverageindex))
            generation = 0
            for i in $chunk
                empty!(got)
                generation += 1
                _checkleafcoverage!(
                    findings,
                    findinglock,
                    got,
                    tree,
                    plan,
                    nearpred,
                    trait,
                    maxfindings,
                    receivingtree_,
                    translatingtree_,
                    planlevels,
                    expected,
                    leafnodes[i],
                    nearlists,
                    coverageindex,
                    stamps,
                    generation,
                )
            end
        end
    end
    return nothing
end

# Coverage partition check ##################################################################
#
# Every receiving leaf must see every translating value exactly once, so the check is
# inherently O(leaves x values). What it does NOT have to be is O(leaves x values x log
# values): sorting each leaf's gathered values and comparing the sorted vector against
# `expected` was 84% of the coverage pass (measured per translating plan at 100k points: 6.7s
# of sort+compare against 0.26s of near-node scanning and 1.0s of far-value gathering).
#
# Stamping replaces that with one linear pass. Each expected value owns a slot; a leaf marks
# the slots it covers with its own generation number, so a slot already carrying that number
# is a duplicate, and a covered count short of `length(expected)` is a gap. No sort, no
# comparison vector, and no per-leaf allocation.
#
# The findings are identical, including which value a duplicate finding
# names. `_firstduplicate` reported the smallest duplicated value because it read a sorted
# vector, so the loop below takes the minimum rather than the first one it happens to meet.

"""
    _CoverageIndex

Slot assignment for the coverage stamp: expected value `v` occupies slot `v - offset`.

`isexpected` is false for ids inside the span that no leaf actually stores, which is how a
covered value that is not an expected value gets reported rather than silently accepted.
"""
struct _CoverageIndex
    offset::Int
    isexpected::Vector{Bool}
    nexpected::Int
end

_nslots(index::_CoverageIndex) = length(index.isexpected)

"""
    _coverageindex(expected)

Build the stamp slots for `expected` (sorted), or `nothing` to keep the sorting path.

`nothing` is returned when the value ids span far more than they populate. The slot array is
proportional to the span, so a sparse one would cost more than the sort it replaces. Also when
`expected` contains a repeat, since then coverage is a multiset question that stamping cannot
answer.
"""
function _coverageindex(expected)
    isempty(expected) && return nothing
    span = last(expected) - first(expected) + 1
    span > 8 * length(expected) + 1024 && return nothing

    offset = first(expected) - 1
    isexpected = zeros(Bool, span)
    for value in expected
        slot = value - offset
        isexpected[slot] && return nothing
        isexpected[slot] = true
    end
    return _CoverageIndex(offset, isexpected, length(expected))
end

"""
    _checkleafcoverage!(findings, findinglock, got, tree, plan, nearpred, trait,
        maxfindings, receivingtree, translatingtree, planlevels, expected, leaf, nearlists,
        coverageindex, stamps, generation)

Check the coverage partition for one receiving leaf using `got` as reusable
scratch storage.
"""
function _checkleafcoverage!(
    findings,
    findinglock,
    got,
    tree,
    plan,
    nearpred,
    trait,
    maxfindings,
    receivingtree,
    translatingtree,
    planlevels,
    expected,
    leaf,
    nearlists,
    coverageindex,
    stamps,
    generation,
)
    _findingsfull(findings, findinglock, maxfindings) && return nothing

    _appendnearvalues!(
        got, tree, receivingtree, translatingtree, leaf, nearpred, trait, nearlists
    )
    for node in Iterators.flatten(((leaf,), ParentUpwardsIterator(receivingtree, leaf)))
        nodelevel = level(receivingtree, node)
        nodelevel in planlevels || continue
        # `plan[node, level]` (not `translatingnodes`): the indexing form returns an empty
        # vector for a node the plan does not address, where the accessor would `KeyError`.
        for translatingnode in plan[node, nodelevel]
            appendvalues!(got, translatingtree, translatingnode)
        end
    end
    duplicate, covered = _tallycoverage!(got, expected, coverageindex, stamps, generation)

    if duplicate !== nothing
        _pushfinding!(
            findings,
            findinglock,
            AdmissibilityFinding(
                :error,
                :coverageduplicate,
                level(receivingtree, leaf),
                leaf,
                missing,
                missing,
                missing,
                missing,
                "leaf $(leaf) covers value $(duplicate) more than once across its near field " *
                "and the far values its translation plan schedules",
            ),
            maxfindings,
        )
    elseif covered != length(expected)
        _pushfinding!(
            findings,
            findinglock,
            AdmissibilityFinding(
                :error,
                :coveragegap,
                level(receivingtree, leaf),
                leaf,
                missing,
                missing,
                missing,
                missing,
                "leaf $(leaf) near field plus plan-scheduled far values do not partition the " *
                "translating side ($(length(got)) covered vs $(length(expected)) expected)",
            ),
            maxfindings,
        )
    end
    return nothing
end

"""
    _appendnearvalues!(out, tree, receivingtree, translatingtree, leaf, nearpred, trait,
                       nearlists)

Append all near-field values for one receiving leaf.

For block trees, values are collected from the opposite side of the block tree.
"""
function _appendnearvalues!(
    out,
    tree,
    receivingtree,
    translatingtree,
    leaf,
    nearpred,
    ::AbstractTreeTrait,
    nearlists,
)
    return appendnearnodevalues!(
        out, receivingtree, leaf; isnear=nearpred, nearlists=nearlists
    )
end

# For a BlockTree the values live on the OPPOSITE side: a test leaf's near set is trial-side
# values. The two-tree `nearnodevalues` takes `(valuetree, referencetree, referencenode)`, so the
# translating tree comes first.
function _appendnearvalues!(
    out, tree, receivingtree, translatingtree, leaf, nearpred, ::isBlockTree, nearlists
)
    if receivingtree === testtree(tree)
        flippednear =
            (trialtree, testtree, trialnode, testnode) ->
                nearpred(testtree, trialtree, testnode, trialnode)
        return appendnearnodevalues!(
            out, translatingtree, receivingtree, leaf; isnear=flippednear
        )
    else
        return appendnearnodevalues!(
            out, translatingtree, receivingtree, leaf; isnear=nearpred
        )
    end
end

"""
    _allvalues(tree)

Return all stored values in the leaves of `tree`, sorted for coverage
comparison.
"""
function _allvalues(tree)
    out = Int[]
    for leaf in leaves(tree)
        appendvalues!(out, tree, leaf)
    end
    sort!(out)
    return out
end

"""
    _tallycoverage!(got, expected, coverageindex, stamps, generation)

Return `(duplicate, covered)` for one leaf's gathered values.

`duplicate` is the smallest value covered more than once, or `nothing`. `covered` is the
number of distinct expected values covered, or `-1` when `got` holds a value that is not an
expected value at all. `-1` never equals `length(expected)`, so the caller reports the same
coverage gap the sorted comparison used to.

Falls back to sorting when `coverageindex` is `nothing`; the two paths agree by construction,
including which duplicate they name.
"""
function _tallycoverage!(got, expected, coverageindex::Nothing, stamps, generation)
    sort!(got)
    # `length(got)` would not do: two vectors of equal length can still differ.
    return _firstduplicate(got), got == expected ? length(expected) : -1
end

function _tallycoverage!(got, expected, coverageindex::_CoverageIndex, stamps, generation)
    duplicate = nothing
    unknown = false
    covered = 0
    for value in got
        slot = value - coverageindex.offset
        if slot < 1 || slot > _nslots(coverageindex) || !coverageindex.isexpected[slot]
            unknown = true
        elseif stamps[slot] == generation
            # Smallest, not first seen: `_firstduplicate` read a sorted vector, and the
            # reported value is part of the finding's message.
            duplicate = isnothing(duplicate) ? value : min(duplicate, value)
        else
            stamps[slot] = generation
            covered += 1
        end
    end
    return duplicate, unknown ? -1 : covered
end

"""
    _firstduplicate(sortedvalues)

Return the first duplicated value in sorted storage, or `nothing` when all
values are unique.
"""
function _firstduplicate(sortedvalues)
    for i in 2:length(sortedvalues)
        sortedvalues[i] == sortedvalues[i - 1] && return sortedvalues[i]
    end
    return nothing
end

"""
    _pushfinding!(findings, findinglock, finding, maxfindings)

Append `finding` unless the report has reached `maxfindings`.

Returns whether the finding was stored.
"""
function _pushfinding!(findings, findinglock, finding, maxfindings)
    lock(findinglock)
    try
        length(findings) >= maxfindings && return false
        push!(findings, finding)
        return true
    finally
        unlock(findinglock)
    end
end

"""
    _findingsfull(findings, findinglock, maxfindings)

Return whether no more findings should be collected.
"""
function _findingsfull(findings, findinglock, maxfindings)
    lock(findinglock)
    try
        return length(findings) >= maxfindings
    finally
        unlock(findinglock)
    end
end

"""
    _raise(report)

Emit warnings and throw on errors according to an [`AdmissibilityReport`](@ref).
"""
function _raise(report::AdmissibilityReport)
    warnings = filter(f -> f.severity === :warning, report.findings)
    isempty(warnings) ||
        @warn "checkadmissibility: $(length(warnings)) marginal far pair(s)" first = first(
            warnings
        )
    if !report.ok
        errors = filter(f -> f.severity === :error, report.findings)
        Base.error(
            "checkadmissibility failed with $(length(errors)) error(s). First: " *
            "[$(first(errors).kind)] $(first(errors).message)",
        )
    end
    return nothing
end
