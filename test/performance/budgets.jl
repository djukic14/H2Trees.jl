# Documented allocation budgets for the performance-contract suite.
#
# `allocation_ratio = allocated_bytes(build_workload) / Base.summarysize(result_tree)`.
#
# The plan-level targets are the aspirational numbers from the performance-contracts refactor
# plan: `<= 10` (hard budget) initially, `<= 5` (stretch) after further optimization.
#
# CAVEAT since the Hilbert level-major renumbering landed: the bulk-built families
# (`TwoNTree`/`BlockTree`/`SimpleHybridTree`) no longer meet the `<=10x` plan budget, and the
# ratio is no longer comparable to those plan numbers *for them*. Nothing got slower; their
# absolute construction allocation is essentially unchanged (measured +1.9%), but the ratio's
# DENOMINATOR shrank ~8x, because a finished tree is now that much smaller. See the per-family
# note below. The remaining families are built by other paths and are still directly comparable.
const ALLOCATION_RATIO_PLAN_HARD_BUDGET = 10.0
const ALLOCATION_RATIO_PLAN_STRETCH = 5.0

# A tree built from 10x more points should not need a much larger allocation-per-final-byte ratio.
# Measured growth (small -> large, `PERF_SMALL_N` -> `PERF_LARGE_N` points) across every family
# below was <= 1.7x; this is a generous ceiling on top of that, not a tight regression trap.
const ALLOCATION_RATIO_SCALING_TOLERANCE = 2.2

# Per-family hard budgets actually asserted by `allocations.jl`, each set with modest headroom
# above the measured baseline (recorded next to each constant; measured on Julia 1.10/1.12,
# `--startup-file=no`, `PERF_SMALL_N=200`/`PERF_LARGE_N=2000` points, 3D unless noted). These are
# real regression guards, not just documentation; they will fail if construction gets
# meaningfully worse than what was measured when they were written.
#
# `TwoNTree`/`BlockTree`/`SimpleHybridTree` allocation history, condensed (full detail, including
# the intermediate fixes and their individual measurements, is in git history and the
# [[project_h2trees_performance_contracts]]/[[project_h2trees_bulkbuild_consolidation]] memories):
#
#   - HISTORICAL PROBLEM: the *default* `TwoNTreeBuilder` (`minvalues=70`) needed each point's
#     stopping level before it could insert anything. This used to mean building a full
#     one-point-per-leaf comparison tree via the same per-point trickle-down insertion the real
#     tree used, THEN separately inserting every point one at a time into the real tree:
#     effectively building two trees to build one, the second of which was itself built the slow
#     way. Original measured ratio: ~115-200x.
#   - CURRENT STATE: `bulkbuildtwontree` (`TwoNTree.jl`) is now the ONLY production construction
#     path. The old point-by-point insertion machinery, the comparison-tree-based stopping-level
#     computation, and the legacy `f(tree, node, value)` protrusion call shape it all depended on
#     have been REMOVED from `src/` entirely, not just superseded. It builds real nodes directly in
#     a single recursive pass, replacing the duplicated derivation above, and its own remaining
#     per-internal-node scratch allocations (a `centers` vector recomputed via `sectorcenter`
#     instead of stored; `nodes`/bucket vectors `sizehint!`ed instead of growing via repeated
#     `push!`-doubling; Hilbert sibling ordering via a stack-allocated `MVector` insertion sort
#     instead of heap-allocating a `Vector` to `sort`) have since been cut too. Measured max:
#     ~5.4x (N=3, large), inside the plan's `<=10x` hard budget and close to the
#     `<=5x` stretch target. See `test/trees/test_bulkbuild.jl`/`test_blocktree_bulkbuild.jl`/
#     `test_uniformseparationdepth.jl` for the correctness coverage this consolidation relies on.
#   - DENOMINATOR CHANGE (Hilbert level-major renumbering): the measured ratio for these three
#     families jumped from ~5.4x to ~40-48x, and this is NOT a construction regression. The
#     numerator barely moved (113,288 -> 115,448 bytes on the 2000-point default 3D build,
#     +1.9%); `summarysize(tree)` fell ~8x (363,584 -> 49,424 bytes). Cause: the bulk builder
#     `sizehint!`s its node vector to a loose upper bound (`min(2*npoints, 1 + 2^N * npoints)`)
#     and the finished tree used to retain that whole oversized backing buffer forever: 4000
#     slots for a 73-node tree on that workload, and 4000 slots for a ONE-node tree on the
#     single-node path. `_hilbertlevelmajornodes` now returns an exactly-sized vector, so the
#     waste is gone. The ratio therefore got numerically worse while actual memory use got ~8x
#     better; the budgets below were re-measured against the new (tight) denominator. Treat the
#     absolute numbers as incomparable to the pre-renumbering ones, but they still guard
#     regressions: a 2x construction-allocation regression would still blow through them.
#   - NUMERATOR CUT (`_uniformseparationdepth` rewrite): construction allocation then dropped
#     ~2.5x (100k 3D points: 68.8MB -> 27.2MB, 36.5ms -> 21.7ms), taking these three families
#     from ~40-48x back to ~10-16x. Two causes, both in the depth-cap scan that used to dominate
#     the build (~65% of all construction bytes): it allocated a fresh `Vector{Vector{Int}}` of
#     per-sector buckets at EVERY node (now one reusable scratch buffer, partitioned in place),
#     and it recursed until every point sat alone in its own cell even though the adaptive path
#     stops a branch at `<= minvalues` (now capped at `minvalues`, which is the depth the builder
#     can actually reach, verified to produce byte-identical trees).
const ALLOCATION_RATIO_BUDGET = Dict(
    # Re-measured after the `_uniformseparationdepth` rewrite (see NUMERATOR CUT above); max
    # across N=1,2,3 x small,large was 15.8x. Headroom kept modest so these stay real regression
    # guards; they would catch the pre-rewrite behaviour returning.
    :TwoNTree => 20.0,           # measured max 15.8x (N=2, large)
    :BlockTree => 20.0,          # measured max 12.4x (small)
    :SimpleHybridTree => 20.0,   # wraps a TwoNTree; same cost profile, measured max 12.3x
    :BoundingBallTree => 11.0,   # measured max 8.8x (large); meets the plan's <=10x
    :KMeansTree => 13.0,         # measured max 10.4x (large); n_threads=1 pinned in the workload.
    # Was 12.6x. The splitwrapper (`ext/H2ParallelKMeansTrees`) built its k-means input matrix
    # via `points[globalpointids]` (a full temporary `Vector{SVector}` immediately fed into
    # `reduce(hcat, ...)` and discarded), read cluster centers back out via
    # `kresult.centers[:, i]` (a temporary column `Vector` per cluster, just to construct an
    # `SVector` from it), and computed radii via `points[partitions[i]]` (another temporary
    # per-partition `Vector{SVector}`, just to iterate it once). None of these three needed to
    # materialize an intermediate array: the k-means matrix is now filled directly by index, the
    # centers are read via `view`, and radii are computed by iterating partition ids directly
    # against `points`; see `_computeradius`/`kmeanswrapper` in
    # `ext/H2ParallelKMeansTrees/H2ParallelKMeansTrees.jl`. Re-measured ~10.6x after a follow-up
    # fix for k-means returning an empty cluster (see that file's comments). The first attempt
    # at that fix reassigned `partitions`/`centers` inside an `if`, which Julia's compiler boxes
    # because both are captured by a closure a few lines below; that silently cost ~24KB/call
    # (ratio back up to ~12.7x) despite the branch never being taken on this workload. Fixed by
    # assigning each exactly once via a ternary instead.
    # allocation-sensitive function grows an `if` that reassigns a variable also used in a
    # comprehension/closure later in the same function.
    :MetisTree => 21.0,          # measured max 17.2x (large); dominated by METIS's own
    # `induced_subgraph` and `BoundingSphere.boundingsphere`'s unavoidable internal `copy(pts)`
    # (its non-mutating entry point always copies before running Welzl's algorithm in place), not
    # by anything H2Trees controls; not revisited here.
    :MetisForest => 12.0,        # measured 9.3x (small only; forest of two small components).
    # `_updatevalues!` (remapping each component's local leaf values back to global point ids
    # after building it as an independent `MetisTree`) used to build a temporary
    # `localtoglobal[values(tree,leaf)]` vector per leaf just to `empty!`/`append!` it straight
    # back in; now remaps every leaf's values vector in place. Real fix, but too small a share of
    # MetisForest's total cost (dominated by the same METIS/BoundingSphere costs as `:MetisTree`)
    # to move the measured ratio outside noise.
)

# Plans build far less than the tree itself (they only ever store node/level index bookkeeping,
# never new geometry), measured ~0.5-1.0x relative to `summarysize(plans)`, already inside the
# stretch target.
const ALLOCATION_RATIO_PLAN_BUDGET = 5.0

# `checkadmissibility` is documented as "a diagnostic, not a hot-path check" (its own docstring),
# so the plan does not set a `<=10x`-style target for it the way it does for construction. Its
# natural "result" (an `AdmissibilityReport`) is tiny regardless of tree size, so the ratio here is
# measured against `summarysize(tree)` instead, a fairer proxy for "work done", since coverage
# checking (`coverage=true`, the default) walks every leaf's near field and the plan's scheduled
# far values. This budget is informational headroom, not a plan-mandated target.
#
# `_checkplancoverage!`'s per-leaf near-field scan (`_checkleafcoverage!`) used to allocate a
# fresh `got::Vector{Int}` per leaf via `nearnodevalues`; it now splits the leaves into up to
# `Threads.nthreads()` contiguous chunks (`_staticchunks`) and spawns one task per chunk
# (`Threads.@spawn`), each reusing its own `appendnearnodevalues!` buffer (see
# `_appendnearvalues!`), emptied between leaves within that chunk, for the task's whole lifetime.
#
# This is deliberately NOT `Threads.@threads :static` with `scratch[Threads.threadid()]`, despite
# that being the more obvious "one buffer per thread" idiom (and what an earlier version of this
# fix used). TWO real problems with that approach were caught under threading and
# realistic call patterns before trusting it, not just by a serial correctness run:
#   - `Threads.threadid()` is not guaranteed to fall within `1:Threads.nthreads()` (Julia's
#     `:default`/`:interactive` thread pools mean the actual bound is `Threads.maxthreadid()`):
#     indexing an `nthreads()`-sized pool by it threw `BoundsError` under `JULIA_NUM_THREADS=4`.
#   - `Threads.@threads :static` itself errors ("cannot be used concurrently or nested") if
#     `checkadmissibility` is ever called from inside the caller's own threaded loop, or if two
#     `checkadmissibility` calls happen to run at once via separate `@spawn`ed tasks; both
#     realistic ways to call a public diagnostic function. `@spawn`ed tasks compose and nest
#     freely, so manual chunking sidesteps this entirely; see
#     `test/plans/test_checkadmissibility.jl`'s "safe to call nested inside, or concurrently
#     with, another threaded region" testset, which pins this.
#
# Measured (`perf_points(3)`, small fixture): Galerkin ~9.1x -> ~3.3x (serial) / ~4.9x (4 threads);
# Petrov/BlockTree ~4.9x -> ~2.0x (serial) / ~2.5x (4 threads).
const ALLOCATION_RATIO_CHECKADMISSIBILITY_BUDGET = 25.0
