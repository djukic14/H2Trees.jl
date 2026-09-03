# Compact PR-comparison benchmark suite for H2Trees.
#
# What is measured here: roughly a dozen representative cases across tree construction, plan
# construction, iterator consumption, interaction-list construction, and one end-to-end workflow.
# `benchmark/ci.jl` compares this suite between a PR's head and base commits; the workflow that
# triggers that comparison is `.github/workflows/Benchmarks.yml`.
#
# Kept intentionally small (roughly a dozen cases, not every method/parameter combination) so a
# base-against-PR comparison run stays fast. See the individual `seconds=`/`samples=` choices below,
# each picked from a warmed-up timing of the underlying operation. Plan construction and
# `checkadmissibility` are inherently much more expensive per call than tree construction (a
# diagnostic walk of the full near/far structure, not a hot-path operation, see their own
# docstrings), so those cases use smaller inputs and fewer samples than tree construction does.
#
# For broader coverage (more sizes, dimensions, distributions, tree families), add an
# `extended.jl` run manually or on a schedule rather than growing this file. See the plan this
# suite was built from for that follow-up.

using BenchmarkTools
using H2Trees
using Random: Xoshiro
using StaticArrays

const SUITE = BenchmarkGroup()

# Deterministic inputs #########################################################################

"""
    detpoints(d, n; seed=0x12345678)

`n` deterministic points in `[0,1]^d`, drawn from a seeded `Xoshiro` RNG so the same points are
generated on every run, which benchmark timings and allocations need to be comparable run to run.
"""
function detpoints(d::Int, n::Int; seed=0x12345678)
    rng = Xoshiro(seed)
    X = rand(rng, d, n)
    return [SVector{d,Float64}(view(X, :, i)) for i in 1:n]
end

# Tree construction #############################################################################
#
# One commonly used builder configuration (the default `TwoNTreeBuilder()`) across 1D/2D/3D, at a
# moderate size and, for 2D/3D, a larger size (>= 50,000 points). Point generation happens here,
# outside the timed expressions below, so each case measures construction only.

const POINTS_1D_10K = detpoints(1, 10_000)
const POINTS_2D_10K = detpoints(2, 10_000)
const POINTS_2D_100K = detpoints(2, 100_000)
const POINTS_3D_10K = detpoints(3, 10_000)
const POINTS_3D_100K = detpoints(3, 100_000)

SUITE["tree construction"]["1D"]["10k"] = @benchmarkable(
    TwoNTree($POINTS_1D_10K), seconds = 1, samples = 20
)
SUITE["tree construction"]["2D"]["10k"] = @benchmarkable(
    TwoNTree($POINTS_2D_10K), seconds = 1, samples = 20
)
SUITE["tree construction"]["2D"]["100k"] = @benchmarkable(
    TwoNTree($POINTS_2D_100K), seconds = 2, samples = 10
)
SUITE["tree construction"]["3D"]["10k"] = @benchmarkable(
    TwoNTree($POINTS_3D_10K), seconds = 1, samples = 20
)
SUITE["tree construction"]["3D"]["100k"] = @benchmarkable(
    TwoNTree($POINTS_3D_100K), seconds = 2, samples = 10
)

# Plan construction ##############################################################################
#
# Trees are built once, outside the timed expressions, so these measure `buildplans` alone. A
# small `minhalfsize`/`minvalues=1` builds a genuinely deep tree; the default `minvalues=70`
# would barely subdivide 2,000 points, making the plan close to trivial.

const PLAN_BUILDER = TwoNTreeBuilder(; minhalfsize=0.02, minvalues=1)
const PLAN_POINTS = detpoints(3, 2_000; seed=0x22345678)
const PLAN_TRIAL_POINTS = detpoints(3, 2_000; seed=0x32345678)
const PLAN_TREE = TwoNTree(PLAN_POINTS; builder=PLAN_BUILDER)
const PLAN_BLOCKTREE = buildtree(
    PLAN_POINTS,
    PLAN_TRIAL_POINTS;
    builder=BlockTreeBuilder(; test=PLAN_BUILDER, trial=PLAN_BUILDER),
)

SUITE["plan construction"]["Galerkin"] = @benchmarkable(
    buildplans($PLAN_TREE; builder=PlanBuilder()), seconds = 3, samples = 10
)
SUITE["plan construction"]["Petrov"] = @benchmarkable(
    buildplans($PLAN_BLOCKTREE; builder=PlanBuilder()), seconds = 3, samples = 10
)

# Iterators #######################################################################################
#
# Full consumption, not merely iterator construction: `foreach` or an explicit sum forces every
# element to actually be produced, so the cost can't be optimized away.

const ITER_TREE = TwoNTree(detpoints(3, 20_000; seed=0x42345678))

function _sumnodes(tree)
    s = 0
    for node in H2Trees.DepthFirstIterator(tree)
        s += node
    end
    return s
end

function _sumnearfar(tree)
    near = 0
    far = 0
    for node in H2Trees.DepthFirstIterator(tree)
        for _ in H2Trees.NearNodeIterator(tree, node)
            near += 1
        end
        for _ in H2Trees.FarNodeIterator(tree, node)
            far += 1
        end
    end
    return near, far
end

SUITE["iterators"]["node traversal"] = @benchmarkable(
    _sumnodes($ITER_TREE), seconds = 1, samples = 20
)
SUITE["iterators"]["leaf traversal"] = @benchmarkable(
    H2Trees.leaves($ITER_TREE), seconds = 1, samples = 20
)
SUITE["iterators"]["near/far interactions"] = @benchmarkable(
    _sumnearfar($ITER_TREE), seconds = 1, samples = 20
)

# Interaction construction #######################################################################
#
# Building near-field value lists and running the admissibility diagnostic, both of which walk the full
# near/far structure the tree and its admissibility criteria establish, so tree shape and the
# near/far predicate directly affect their cost.

function _gathernearfield(tree)
    out = Int[]
    for leaf in H2Trees.leaves(tree)
        appendnearnodevalues!(out, tree, leaf)
    end
    return out
end

const INTERACTION_PLANS = buildplans(PLAN_TREE; builder=PlanBuilder())

SUITE["interaction construction"]["near-field lists"] = @benchmarkable(
    _gathernearfield($ITER_TREE), seconds = 1, samples = 20
)
SUITE["interaction construction"]["checkadmissibility"] = @benchmarkable(
    checkadmissibility($PLAN_TREE, $INTERACTION_PLANS; throw=false),
    seconds = 3,
    samples = 10
)

# Representative end-to-end workflow #############################################################
#
# input -> tree -> plan -> consume interactions, in one timed expression, which catches regressions
# that only show up when the whole pipeline runs together, not in any isolated microbenchmark.

const WORKFLOW_POINTS = detpoints(3, 2_000; seed=0x62345678)

function _endtoendworkflow(points, builder)
    tree = TwoNTree(points; builder=builder)
    plans = buildplans(tree; builder=PlanBuilder())
    return checkadmissibility(tree, plans; throw=false)
end

SUITE["representative workflows"]["input to interactions (3D)"] = @benchmarkable(
    _endtoendworkflow($WORKFLOW_POINTS, $PLAN_BUILDER), seconds = 3, samples = 10
)

# Geometric translation construction ############################################################

const TRANSLATION_POINTS = detpoints(3, 20_000; seed=0x42345678)
const TRANSLATION_TREE = TwoNTree(TRANSLATION_POINTS; builder=PLAN_BUILDER)
const TRANSLATION_PLANS = buildplans(TRANSLATION_TREE; builder=PlanBuilder())
const TRANSLATION_PLAN = H2Trees.translatingplan(
    TRANSLATION_PLANS.trialaggregationplan, TRANSLATION_PLANS.testdisaggregationplan
)

# Trees and plans are built above, outside the timed expressions, so each case measures
# `translations` alone, which is what the comparison between them is about.
SUITE["translation construction"]["direction only"] = @benchmarkable(
    H2Trees.translations(
        $TRANSLATION_TREE, $TRANSLATION_PLAN, H2Trees.DirectionInvariancePerLevel()
    ),
    seconds = 3,
    samples = 5
)
SUITE["translation construction"]["symmetry: opposite"] = @benchmarkable(
    H2Trees.translations(
        $TRANSLATION_TREE,
        $TRANSLATION_PLAN,
        H2Trees.SymmetryDirectionInvariancePerLevel(H2Trees.OppositeSymmetry()),
    ),
    seconds = 3,
    samples = 5
)
SUITE["translation construction"]["symmetry: axis-preserving"] = @benchmarkable(
    H2Trees.translations(
        $TRANSLATION_TREE,
        $TRANSLATION_PLAN,
        H2Trees.SymmetryDirectionInvariancePerLevel(H2Trees.AxisPreservingSymmetry(3)),
    ),
    seconds = 3,
    samples = 5
)
SUITE["translation construction"]["symmetry: full lattice"] = @benchmarkable(
    H2Trees.translations(
        $TRANSLATION_TREE,
        $TRANSLATION_PLAN,
        H2Trees.SymmetryDirectionInvariancePerLevel(H2Trees.FullLatticeSymmetry()),
    ),
    seconds = 3,
    samples = 5
)

# Symmetry reduction when the two block-tree roots are lattice-aligned but not coincident.

const BLOCK_TRANSLATION_POINTS = detpoints(3, 2_000; seed=0x52345678)
const BLOCK_TRANSLATION_BUILDER = TwoNTreeBuilder(; minhalfsize=0.02, minvalues=1)

# One ROOT halfsize, read off a tree built from the same points, so the shift below is a lattice
# vector at every level rather than only at the finest.
const BLOCK_ROOT_HALFSIZE = let
    tree = TwoNTree(BLOCK_TRANSLATION_POINTS; builder=BLOCK_TRANSLATION_BUILDER)
    H2Trees.halfsize(tree, H2Trees.root(tree))
end

function blocktranslationtree(shift)
    return buildtree(
        [p + SVector{3,Float64}(shift, 0.0, 0.0) for p in BLOCK_TRANSLATION_POINTS],
        BLOCK_TRANSLATION_POINTS;
        builder=BlockTreeBuilder(;
            test=BLOCK_TRANSLATION_BUILDER, trial=BLOCK_TRANSLATION_BUILDER
        ),
    )
end

function blocktranslationplan(tree)
    plans = buildplans(tree; builder=PlanBuilder())
    return H2Trees.translatingplan(plans.trialaggregationplan, plans.testdisaggregationplan)
end

const BLOCK_COINCIDENT_TREE = blocktranslationtree(0.0)
const BLOCK_COINCIDENT_PLAN = blocktranslationplan(BLOCK_COINCIDENT_TREE)
const BLOCK_OFFSET_TREE = blocktranslationtree(BLOCK_ROOT_HALFSIZE)
const BLOCK_OFFSET_PLAN = blocktranslationplan(BLOCK_OFFSET_TREE)

# Trees and plans are built outside the timed expressions, as in the single-tree group, so each case
# measures `translations` alone.
for (rootlabel, tree, plan) in (
    ("coincident", BLOCK_COINCIDENT_TREE, BLOCK_COINCIDENT_PLAN),
    ("lattice offset", BLOCK_OFFSET_TREE, BLOCK_OFFSET_PLAN),
)
    SUITE["block translation construction"][rootlabel]["direction only"] = @benchmarkable(
        H2Trees.translations($tree, $plan, H2Trees.DirectionInvariancePerLevel()),
        seconds = 3,
        samples = 5
    )
    SUITE["block translation construction"][rootlabel]["symmetry: opposite"] = @benchmarkable(
        H2Trees.translations(
            $tree,
            $plan,
            H2Trees.SymmetryDirectionInvariancePerLevel(H2Trees.OppositeSymmetry()),
        ),
        seconds = 3,
        samples = 5
    )
    SUITE["block translation construction"][rootlabel]["symmetry: axis-preserving"] = @benchmarkable(
        H2Trees.translations(
            $tree,
            $plan,
            H2Trees.SymmetryDirectionInvariancePerLevel(H2Trees.AxisPreservingSymmetry(3)),
        ),
        seconds = 3,
        samples = 5
    )
end
