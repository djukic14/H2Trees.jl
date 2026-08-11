# Locality quality measurement for the Hilbert node ordering.
#
# Manual locality report, not a timing benchmark or correctness test.
#
#     julia --project=benchmark benchmark/locality.jl
#
# WHAT IT MEASURES
#
# The motivation for ordering nodes along a Hilbert curve is that a fixed-size chunk of leaf ids
# should cover a compact region, so the near-field data its leaves need is fetched once and
# reused. For a chunk `C` of leaves that is:
#
#     |union of near(B) for B in C| / |C|
#
# Distinct near boxes touched per leaf. Lower means better near-field reuse.
#
# ORDERS COMPARED
#
#   - `dfs`    : `leaves(tree)` as returned.
#   - `byid`   : `sort(leaves(tree))`, level-major then Hilbert within each level.
#   - `random` : a shuffle, as the no-locality control.
#
# WHAT THE NUMBERS SHOWED WHEN THIS WAS WRITTEN
#
#   - Balanced trees: `dfs` and `byid` coincide, and both beat `random` by ~1.6-4.3x.
#   - Adaptive trees (leaves on several levels): `byid` beats `dfs` by ~1.14-1.63x. Near-field
#     coupling is predominantly same-level, so grouping leaves by level puts same-scale boxes
#     with heavily overlapping near fields in one chunk, whereas depth-first interleaves scales.
#     This is where the renumbering actually pays off.
#
# Re-run this after any change to node ordering, near/far predicates, or the bulk builder's
# partitioning; a large regression in the `byid` column means chunked near-field work got worse.

using H2Trees
using Random: Xoshiro, shuffle
using StaticArrays
using Printf

"""
    chunkcost(tree, order, chunksize)

Mean over chunks of `|union of near(B) for B in C| / |C|`. Lower is better.
"""
function chunkcost(tree, order, chunksize)
    total = 0.0
    nchunks = 0
    seen = Set{Int}()
    for chunk in Iterators.partition(order, chunksize)
        empty!(seen)
        for leaf in chunk
            for nearnode in H2Trees.NearNodeIterator(tree, leaf)
                push!(seen, nearnode)
            end
        end
        total += length(seen) / length(chunk)
        nchunks += 1
    end
    return total / nchunks
end

# Deterministic geometries. The balanced group isolates "does the ordering help at all"; the
# adaptive group is where the level-major layout differs from depth-first at all.
function geometries()
    cases = Tuple{String,Symbol,Vector{<:SVector},Any}[]

    side = 8
    push!(
        cases,
        (
            "cartesian 8^3",
            :balanced,
            [
                SVector(i + 0.5, j + 0.5, k + 0.5) for i in 0:(side - 1) for
                j in 0:(side - 1) for k in 0:(side - 1)
            ],
            (; minhalfsize=0.5, minvalues=0),
        ),
    )

    rng = Xoshiro(5)
    sphere = SVector{3,Float64}[]
    for _ in 1:4000
        v = randn(rng, 3)
        v ./= sqrt(sum(abs2, v))
        push!(sphere, SVector(v...))
    end
    push!(cases, ("sphere surface", :balanced, sphere, (; minhalfsize=0.08, minvalues=0)))

    rng = Xoshiro(11)
    twoscale = SVector{3,Float64}[]
    for _ in 1:3000
        push!(twoscale, SVector((0.03 .* randn(rng, 3))...))
    end
    for _ in 1:600
        push!(twoscale, SVector(randn(rng, 3)...))
    end
    push!(cases, ("two-scale core+halo", :adaptive, twoscale, (; minvalues=32)))

    rng = Xoshiro(12)
    clusters = SVector{3,Float64}[]
    for c in ((0, 0, 0), (3, 0, 0), (0, 3, 0), (3, 3, 3))
        for _ in 1:(c == (0, 0, 0) ? 2500 : 300)
            push!(clusters, SVector((c .+ 0.15 .* randn(rng, 3))...))
        end
    end
    push!(cases, ("uneven clusters", :adaptive, clusters, (; minvalues=32)))

    rng = Xoshiro(13)
    mixed = SVector{3,Float64}[]
    for _ in 1:2000
        v = randn(rng, 3)
        v ./= sqrt(sum(abs2, v))
        push!(mixed, SVector(v...))
    end
    for _ in 1:1500
        push!(mixed, SVector((0.05 .* randn(rng, 3))...))
    end
    push!(cases, ("surface + volume", :adaptive, mixed, (; minvalues=32)))

    return cases
end

function main()
    @printf(
        "%-22s %-6s %-8s %-8s %-8s %-10s\n",
        "geometry",
        "chunk",
        "dfs",
        "byid",
        "random",
        "byid vs dfs"
    )
    println(repeat("-", 70))

    for (name, kind, points, builderkw) in geometries()
        tree = buildtree(points; builder=TwoNTreeBuilder(; builderkw...))
        leaves = H2Trees.leaves(tree)
        byid = sort(copy(leaves))
        random = shuffle(Xoshiro(7), copy(leaves))

        leaflevels = sort(unique(H2Trees.level(tree, leaf) for leaf in leaves))
        contiguous = leaves == collect(minimum(leaves):maximum(leaves))
        @printf(
            "%s [%s]: %d leaves, leaf levels %s, contiguous ids: %s\n",
            name,
            kind,
            length(leaves),
            leaflevels,
            contiguous
        )

        for chunksize in (32, 64, 128)
            chunksize > length(leaves) && continue
            d = chunkcost(tree, leaves, chunksize)
            b = chunkcost(tree, byid, chunksize)
            r = chunkcost(tree, random, chunksize)
            @printf(
                "%-22s %-6d %-8.3f %-8.3f %-8.3f %-10s\n",
                "",
                chunksize,
                d,
                b,
                r,
                @sprintf("%.2fx", d / b)
            )
        end
        println()
    end
    return nothing
end

main()
