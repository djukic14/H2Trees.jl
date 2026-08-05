# Small wrapped functions exercising construction, iterator consumption, plan construction, and
# diagnostics -- each one is what `inference.jl`/`dispatch.jl`/`allocations.jl` actually calls, so
# a workload here is exactly the unit that gets a budget/inference/dispatch check.

using H2Trees
using StaticArrays
using Random
using BoundingSphere: boundingsphere
using Metis
using ParallelKMeans

# Tree construction ########################################################################

function perf_buildtwontree(points; minvalues::Int=70)
    return buildtree(points; builder=TwoNTreeBuilder(; minvalues=minvalues))
end

function perf_buildblocktree(testpoints, trialpoints; minvalues::Int=70)
    return buildtree(
        testpoints,
        trialpoints;
        builder=BlockTreeBuilder(;
            test=TwoNTreeBuilder(; minvalues=minvalues),
            trial=TwoNTreeBuilder(; minvalues=minvalues),
        ),
    )
end

# A low `minvalues` paired with an explicit nonzero `minhalfsize` builds genuinely deep,
# well-separated trees for the plan/admissibility fixtures, where the default `minvalues=70`
# would (on a few hundred points) barely subdivide at all.
function perf_buildplantree(points; minhalfsize::Float64=0.05)
    return buildtree(
        points; builder=TwoNTreeBuilder(; minhalfsize=minhalfsize, minvalues=1)
    )
end

function perf_buildplanblocktree(testpoints, trialpoints; minhalfsize::Float64=0.05)
    return buildtree(
        testpoints,
        trialpoints;
        builder=BlockTreeBuilder(;
            test=TwoNTreeBuilder(; minhalfsize=minhalfsize, minvalues=1),
            trial=TwoNTreeBuilder(; minhalfsize=minhalfsize, minvalues=1),
        ),
    )
end

function _spreadaxis(points, ids)
    N = length(points[first(ids)])
    bestaxis, bestspread = 1, -Inf
    for d in 1:N
        lo, hi = extrema(points[i][d] for i in ids)
        spread = hi - lo
        if spread > bestspread
            bestaxis, bestspread = d, spread
        end
    end
    return bestaxis
end

# Deterministic median-chunk splitter: no RNG, splits `globalpointids` into `numsplits`
# contiguous groups after sorting along whichever axis currently has the greatest spread.
function perf_chunksplitter(points, globalpointids, level, numsplits; kwargs...)
    axis = _spreadaxis(points, globalpointids)
    sortedids = sort(globalpointids; by=i -> points[i][axis])
    n = length(sortedids)
    base, rem = divrem(n, numsplits)
    partitions = Vector{Int}[]
    idx = 1
    for c in 1:numsplits
        len = base + (c <= rem ? 1 : 0)
        len == 0 && continue
        push!(partitions, sortedids[idx:(idx + len - 1)])
        idx += len
    end
    N = length(points[first(globalpointids)])
    T = eltype(points[first(globalpointids)])
    centers = Vector{SVector{N,T}}(undef, length(partitions))
    radii = Vector{T}(undef, length(partitions))
    for i in eachindex(partitions)
        c, r = boundingsphere(points[partitions[i]])
        centers[i] = SVector(c...)
        radii[i] = r
    end
    return partitions, centers, radii
end

function perf_buildboundingballtree(points; numsplits::Int=4, minvalues::Int=numsplits)
    return buildtree(
        points;
        builder=BoundingBallTreeBuilder(;
            splitter=perf_chunksplitter, numsplits=numsplits, minvalues=minvalues
        ),
    )
end

function perf_buildkmeanstree(points; numberofclusters::Int=4)
    # `n_threads=1` pins ParallelKMeans to a single thread regardless of the ambient
    # `JULIA_NUM_THREADS`: its internal parallelism allocates thread-local buffers per task, so the
    # measured allocation_ratio would otherwise depend on how many threads happen to be available
    # rather than on construction itself -- exactly the kind of environment-dependence this suite
    # is meant to avoid (see `KMeansTreeBuilder`'s own seeded-`rng` default for the same reasoning
    # applied to determinism).
    return buildtree(
        points;
        builder=KMeansTreeBuilder(;
            numberofclusters=numberofclusters,
            splitterkwargs=(; n_threads=1, rng=Random.MersenneTwister(1234)),
        ),
    )
end

function perf_buildmetistree(points, graph, weights; numdivisions::Int=4)
    return buildtree(
        points, graph, weights; builder=MetisTreeBuilder(; numdivisions=numdivisions)
    )
end

function perf_buildmetisforest(points, graph, weights; numdivisions::Int=4)
    return buildforest(
        points,
        graph,
        weights;
        builder=MetisForestBuilder(;
            treebuilder=MetisTreeBuilder(; numdivisions=numdivisions)
        ),
    )
end

function perf_buildsimplehybridtree(points; minvalues::Int=70)
    tree = perf_buildtwontree(points; minvalues=minvalues)
    # The golden-ratio point cloud can isolate a point into its own leaf as shallow as level 2
    # regardless of point count (a sparse octant separates fast), so any hybrid level below the
    # shallowest leaf is fixture-dependent and fragile. Doubling the root's own halfsize is always
    # above every level's halfsize, so `SimpleHybridTree` falls back to `hybridlevel = root(tree)`
    # -- the one choice guaranteed valid for any tree shape, since the root is never a leaf in a
    # tree with more than one node.
    hybridhalfsize = 2 * maximum(H2Trees.halfsizes(tree))
    return buildtree(tree; builder=SimpleHybridTreeBuilder(; hybridhalfsize=hybridhalfsize))
end

# Iterator consumption ######################################################################

function perf_depthfirstsum(tree)
    s = 0
    for node in H2Trees.DepthFirstIterator(tree)
        s += node
    end
    return s
end

function perf_childrensum(tree)
    s = 0
    for node in H2Trees.DepthFirstIterator(tree)
        for child in H2Trees.ChildIterator(tree, node)
            s += child
        end
    end
    return s
end

function perf_parentupwardssum(tree)
    s = 0
    for leaf in H2Trees.leaves(tree)
        for p in H2Trees.ParentUpwardsIterator(tree, leaf)
            s += p
        end
    end
    return s
end

function perf_levelsum(tree)
    s = 0
    for level in H2Trees.levels(tree)
        for node in H2Trees.LevelIterator(tree, level)
            s += node
        end
    end
    return s
end

function perf_nearfarsum(tree)
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

# Plans ######################################################################################

function perf_buildgalerkinplans(tree)
    return buildplans(tree; builder=PlanBuilder())
end

function perf_buildpetrovplans(blocktree)
    return buildplans(blocktree; builder=PlanBuilder())
end

function perf_checkadmissibility(tree, plans)
    return checkadmissibility(tree, plans; throw=false)
end
