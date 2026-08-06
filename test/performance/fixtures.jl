# Deterministic fixtures for the performance-contract suite.
#
# Everything here must be reproducible byte-for-byte across runs (no `Random.GLOBAL_RNG`, no
# wall-clock seeds) so allocation/inference/dispatch measurements are comparable run to run.

using StaticArrays
using Graphs

"""
    detpoints(N, n; scale=1.0, offset=zero(SVector{N,Float64}))

`n` deterministic points in `N` dimensions, generated from an irrational-rotation lattice
(golden-ratio equidistribution) so coordinates never land exactly on a shared axis value or
coincide -- both of which would trigger degenerate tree-construction edge cases unrelated to
what these tests measure.
"""
function detpoints(
    N::Int, n::Int; scale::Float64=1.0, offset::SVector=SVector(ntuple(_ -> 0.0, N))
)
    golden = (sqrt(5.0) - 1.0) / 2.0
    pts = Vector{SVector{N,Float64}}(undef, n)
    for i in 1:n
        pts[i] = SVector(ntuple(d -> scale * (mod(i * golden^d, 1.0) - 0.5), N)) + offset
    end
    return pts
end

const PERF_SMALL_N = 200
const PERF_LARGE_N = 2000

perf_points(N::Int; large::Bool=false) = detpoints(N, large ? PERF_LARGE_N : PERF_SMALL_N)
function perf_trial_points(N::Int; large::Bool=false)
    return detpoints(
        N, large ? PERF_LARGE_N : PERF_SMALL_N; offset=SVector(ntuple(_ -> 3.0, N))
    )
end

"""
    perf_graph(n)

An `n`-vertex connected graph (a cycle) for `MetisTree` fixtures -- one vertex per point,
guaranteed connected so a single `MetisTree` (rather than a multi-component forest) is built.
"""
perf_graph(n::Int) = Graphs.cycle_graph(n)

"""
    perf_forest_graph(n1, n2)

Two disjoint connected components (cycles of size `n1` and `n2`, no edges between them) for
`MetisForest` fixtures, which require more than one connected component to be meaningful.
"""
function perf_forest_graph(n1::Int, n2::Int)
    return Graphs.blockdiag(Graphs.cycle_graph(n1), Graphs.cycle_graph(n2))
end

perf_weights(n::Int) = ones(Int, n)
