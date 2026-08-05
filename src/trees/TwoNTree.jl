"""
    TwoNTree{N,D,T} <: H2ClusterTree

Cluster tree whose nodes are bounded by axis-aligned boxes.

`N` is the ambient dimension, `D` is the node-data type, and `T` is the
coordinate/halfsize type. The cached [`TreeIndex`](@ref) is stored in a `Ref` so
`rebuildtreeindex!` can replace it after topology-changing operations.
"""
struct TwoNTree{N,D,T,I} <: H2ClusterTree
    nodes::Vector{Node{D}}
    root::Int
    center::SVector{N,T}
    halfsize::T
    index::I
end

"""
    TwoNTree(nodes, root, center, halfsize, nodesatlevel)

Construct a `TwoNTree` from already materialized nodes and level storage, then
rebuild the cached tree index.
"""
function TwoNTree(
    nodes::Vector{Node{D}},
    root::Int,
    center::SVector{N,T},
    halfsize::T,
    nodesatlevel::Vector{Vector{Int}},
) where {N,D,T}
    minlevel = isempty(nodesatlevel) ? 0 : level(first(nodes).data)
    maxlevel = isempty(nodesatlevel) ? 0 : minlevel + length(nodesatlevel) - 1
    index = Ref(TreeIndex(nodesatlevel, Int[], Int[], minlevel, maxlevel))
    tree = TwoNTree{N,D,T,Base.RefValue{TreeIndex}}(nodes, root, center, halfsize, index)
    return rebuildtreeindex!(tree)
end

function TwoNTree{N,D,T}(
    nodes::Vector{Node{D}},
    root::Int,
    center::SVector{N,T},
    halfsize::T,
    nodesatlevel::Vector{Vector{Int}},
) where {N,D,T}
    return TwoNTree(nodes, root, center, halfsize, nodesatlevel)
end

"""
    tree(node)

Return the node object stored at global node id `node`.
"""
function (tree::H2ClusterTree)(node::Int)
    return tree.nodes[node - H2Trees.root(tree) + 1]
end

"""
    TwoNTree(positions; builder::TwoNTreeBuilder=TwoNTreeBuilder())

Construct a `TwoNTree` from positions.

Prefer the canonical [`buildtree`](@ref) entry point; this method is the
builder-backed constructor it forwards to.
"""
function TwoNTree(positions; builder::TwoNTreeBuilder=TwoNTreeBuilder())
    builder = _resolve_builder_protrusion(builder, positions)
    rootcenter, rootsize = boundingbox(positions)
    # `minhalfsize` may arrive as an untyped default (e.g. the `0` from `TwoNTreeBuilder()`);
    # promote it to the coordinate type so it matches `roothalfsize`/the inner constructor.
    minhalfsize = oftype(rootsize, builder.minhalfsize)
    minlevel = _resolve_minlevel(builder.minlevel, 1)

    # The inner constructor below (`bulkbuildtwontree`) already rebuilds the tree index once as
    # part of construction: rebuilding it again here would just redo that same O(n log n)
    # traversal against an unchanged tree.
    return TwoNTree(
        SVector(rootcenter...),
        positions,
        roothalfsize(rootsize, minhalfsize),
        minhalfsize;
        minlevel=minlevel,
        root=builder.root,
        minvalues=builder.minvalues,
        protrusion=builder.protrusion,
    )
end

"""
    TwoNTree(center, halfsize; minlevel=1, root=1, minvalues=0)

Construct a one-node `TwoNTree`.
"""
function TwoNTree(
    center::SVector{N,T}, halfsize::T; minlevel::Int=1, root::Int=1, minvalues=0
) where {N,T}
    rootnode = Node(BoxData(0, Int[], center, halfsize, minlevel), 0, 0, 0)
    tree = TwoNTree(
        [rootnode],
        root,
        center,
        halfsize,
        Ref(TreeIndex([[root]], [root], [root], minlevel, minlevel)),
    )
    return tree
end

"""
    TwoNTree(center, points, halfsize, minhalfsize; kwargs...)

Build a `TwoNTree` directly from root geometry and point positions.

Most callers should use [`buildtree`](@ref) or `TwoNTree(positions; builder=...)`
instead; this method is the lower-level entry to [`bulkbuildtwontree`](@ref).
"""
function TwoNTree(
    center::SVector{N,T},
    points::AbstractVector{SVector{N,T}},
    halfsize::T,
    minhalfsize::T;
    minlevel::Int=1,
    root::Int=1,
    minvalues=0,
    protrusion=NoProtrusionCheck(),
) where {N,T}
    return bulkbuildtwontree(
        points, center, halfsize, minhalfsize, minlevel, root, minvalues, protrusion
    )
end

"""
    sectorcentersize(point, center, halfsize)

Return `(sector, childcenter, childhalfsize)` for `point` in the box centered at
`center` with halfsize `halfsize`.
"""
function sectorcentersize(pt, ct, hs)
    hs = hs / 2
    bl = pt .> ct
    ct = ifelse.(bl, ct .+ hs, ct .- hs)
    sc = sum(b ? 2^(i - 1) : 0 for (i, b) in enumerate(bl))
    return sc, ct, hs
end

"""
    sectorcenter(sector, center, childhalfsize)

Return the child center for `sector` under a parent centered at `center`.
"""
function sectorcenter(sector::Int, center::SVector{N,T}, childhalfsize::T) where {N,T}
    return SVector(ntuple(i -> if (((sector >> (i - 1)) & 1) == 1)
        center[i] + childhalfsize
    else
        center[i] - childhalfsize
    end, N))
end

"""
    parentcenterminuschildcenter(tree::TwoNTree, child)

Return `center(parent(tree, child)) - center(tree, child)`.
"""
function parentcenterminuschildcenter(tree::TwoNTree{N,D,T}, child::Int) where {N,D,T}
    return parentcenterminuschildcenter(N, sector(tree, child), halfsize(tree, child))
end

function parentcenterminuschildcenter(N, sector, halfsize)
    ds = digits(sector; base=2, pad=N)
    for i in eachindex(ds)
        ds[i] == 0 && (ds[i] = -1)
    end

    return SVector{N}(ds .* -halfsize)
end

"""
    oppositesector(tree::TwoNTree, node)
    oppositesector(N, sector)

Return the sector opposite to `node`'s sector, or to `sector` in `N`
dimensions.
"""
function oppositesector(tree::TwoNTree{N,D,T}, node::Int) where {N,D,T}
    return oppositesector(N, sector(tree, node))
end

function oppositesector(N::Int, sector)
    return (2^N - 1) - sector
end

"""
    _uniformseparationdepth(points, indices, center, halfsize)

Number of recursive octree-bisection levels needed so that every point in `points[indices]` lands
in its own sector, applying the exact same box-sector rule `TwoNTree` itself uses
(`sectorcentersize`) directly to point indices, without building any tree. Returns 0 once
`indices` holds at most one point.

This is a standalone geometric computation on purpose: it must NOT go through `bulkbuildtwontree`
itself, which calls this function directly (as a safety depth cap for its adaptive path) to decide
where to stop: calling back into tree construction here would recurse forever.

Throws `ArgumentError` if `indices` cannot be separated, i.e. two or more of `points[indices]` are
coincident (or differ by less than `Float64` can represent at this scale). `sectorcentersize`'s
sector choice depends only on `center` (not `halfsize`), so coincident points fall in the same
sector at every level no matter how far bisection recurses; `center` converges toward them and
`halfsize` halves toward `0.0`, but the points never separate. Recursing past `halfsize == 0.0`
would never terminate (the geometry stops changing entirely), so that is the signal checked here.
"""
function _uniformseparationdepth(points, indices, center, halfsize)
    length(indices) <= 1 && return 0
    if iszero(halfsize)
        return throw(
            ArgumentError(
                "cannot separate $(length(indices)) coincident points (or points closer than " *
                "floating-point precision can resolve) into distinct boxes; deduplicate the " *
                "input points before building a tree",
            ),
        )
    end

    # We know indices has >1 point, so this call itself performs a real split -- that split
    # counts as one level even when it already fully separates every point into singleton
    # buckets (the recursive calls below then correctly contribute 0 for those).
    #
    # `sc` (from `sectorcentersize`) always lands in `0:2^N-1`, so a plain `Vector` indexed by
    # `sc + 1` covers every sector without the hashing/resizing overhead of a `Dict` -- this
    # function recurses once per tree level, so that overhead is paid on every level of every
    # `bulkbuildtwontree` adaptive-path build.
    nsectors = 1 << length(center)
    buckets = Vector{Vector{Int}}(undef, nsectors)
    childgeometry = Vector{Tuple{typeof(center),typeof(halfsize)}}(undef, nsectors)
    for i in indices
        sc, childcenter, childhalfsize = sectorcentersize(points[i], center, halfsize)
        idx = sc + 1
        if !isassigned(buckets, idx)
            buckets[idx] = Int[]
            childgeometry[idx] = (childcenter, childhalfsize)
        end
        push!(buckets[idx], i)
    end

    maxchilddepth = 0
    for idx in 1:nsectors
        isassigned(buckets, idx) || continue
        childcenter, childhalfsize = childgeometry[idx]
        maxchilddepth = max(
            maxchilddepth,
            _uniformseparationdepth(points, buckets[idx], childcenter, childhalfsize),
        )
    end
    return 1 + maxchilddepth
end

"""
    bulkbuildtwontree(points, rootcenter, roothalfsize, minhalfsize, minlevel, root, minvalues, protrusion)

Build a complete `TwoNTree` in a single recursive pass over grouped point ids.

For each cell, decides whether to stop (create a leaf holding every point index in the cell) or
continue (bucket by [`sectorcentersize`](@ref) and recurse into each non-empty child) using two
different rules depending on `minvalues`/`protrusion`:

  - the general adaptive rule: `length(ids) <= minvalues`, then whether splitting would make any
    child protrude, with a root-only special case (if ANY child of the root would protrude, the
    WHOLE tree stays at the root, not just the offending branch) -- plus two purely geometric
    caps: `halfsize <= minhalfsize`, and the uniform-bisection depth
    [`_uniformseparationdepth`](@ref) computes (also the source of the coincident-point
    `ArgumentError` in this regime -- deduplicate points before building an adaptive tree);
  - `minvalues == 0` with `NoProtrusionCheck`: a uniform-depth rule bounded ONLY by
    `halfsize <= minhalfsize`, with no conformance computation and no coincident-point
    resolution at all (points landing in the same cell all the way down to `minhalfsize` simply
    share a leaf) -- implemented by passing `_bulkbuildnode!`/`_bulkbuildchildren!` a sentinel
    `nlevels = typemax(Int)` so the depth-based stopping condition never fires, skipping
    `_uniformseparationdepth` entirely (see the `nlevels` computation below).

Either way, children are inserted in Hilbert-curve sibling order (`hilbert_positions`/
`hilbert_states`), so `firstchild`/`nextsibling` pointers are deterministic regardless of point
processing order.

`protrusion.compute` is called as `f(center, halfsize, value)`; tree/node-shaped
construction-time protrusion functors are not supported by the bulk builder.
"""
function bulkbuildtwontree(
    points::AbstractVector{SVector{N,T}},
    rootcenter::SVector{N,T},
    roothalfsize::T,
    minhalfsize::T,
    minlevel::Int,
    root::Int,
    minvalues::Int,
    protrusion,
) where {N,T}
    _validateprotrusion(protrusion)
    # `minhalfsize` is the ONLY thing bounding recursion in the `minvalues=0`+`NoProtrusionCheck`
    # uniform-bulk regime below (it deliberately skips `_uniformseparationdepth`'s depth cap) --
    # a negative or non-finite value would mean `halfsize <= minhalfsize` never becomes true,
    # recursing until `halfsize` underflows to exactly `0.0` (or forever, for `NaN`, since even
    # that comparison is always `false`). The public `TwoNTree(positions; builder)` entry point
    # already prevents this indirectly (`roothalfsize` requires a valid `minhalfsize`), but this
    # function is itself callable directly, so it validates on its own rather than trusting every
    # caller to have done so.
    (isfinite(minhalfsize) && minhalfsize >= zero(T)) || throw(
        ArgumentError("minhalfsize must be finite and non-negative, got $minhalfsize")
    )

    maxprotrusion = protrusion isa ProtrusionCheck ? protrusion.max : NaN
    computeprotrusion =
        protrusion isa ProtrusionCheck ? protrusion.compute : ComputeProtrusionFunctor()

    # minvalues=0 + NoProtrusionCheck: bounded purely by minhalfsize, not by whether points can be
    # separated into distinct leaves -- a sentinel nlevels the depth check in
    # `_bulkbuildnode!`/`_bulkbuildchildren!` never reaches, so `_uniformseparationdepth` (which
    # throws for coincident points) is never called at all for this regime. See the docstring above.
    nlevels = if iszero(minvalues) && protrusion isa NoProtrusionCheck
        typemax(Int)
    else
        max(1, _uniformseparationdepth(points, eachindex(points), rootcenter, roothalfsize))
    end

    nodes = Vector{Node{BoxData{N,T}}}()
    # Conservative upper bound on final node count: a full `2^N`-ary tree over `length(points)`
    # leaves has at most `(2^N * length(points) - 1) / (2^N - 1)` nodes; the simpler
    # `(2^N) * length(points)` bound used here is looser but cheap to compute and avoids
    # `sizehint!` itself becoming a source of repeated reallocation on the (much more common)
    # shallower trees.
    sizehint!(nodes, min(2 * length(points), 1 + (1 << N) * length(points)))
    rootoffset = root - 1

    if roothalfsize <= minhalfsize
        push!(
            nodes,
            Node(
                BoxData(0, collect(eachindex(points)), rootcenter, roothalfsize, minlevel),
                0,
                0,
                0,
            ),
        )
        return _finishbulktwontree(nodes, root, rootcenter, roothalfsize)
    end

    nsectors = 1 << N
    childhalfsize = roothalfsize / 2
    buckets = _bucketbysector(points, eachindex(points), rootcenter, roothalfsize, nsectors)

    # Root's own conformance: if ANY child would protrude, the whole tree stays unsplit at the
    # root -- see the docstring above.
    rootconforms = true
    if !isnan(maxprotrusion)
        for idx in 1:nsectors
            isassigned(buckets, idx) || continue
            if _anyprotrudes(
                points,
                buckets[idx],
                sectorcenter(idx - 1, rootcenter, childhalfsize),
                childhalfsize,
                maxprotrusion,
                computeprotrusion,
            )
                rootconforms = false
                break
            end
        end
    end

    if !rootconforms
        push!(
            nodes,
            Node(
                BoxData(0, collect(eachindex(points)), rootcenter, roothalfsize, minlevel),
                0,
                0,
                0,
            ),
        )
        return _finishbulktwontree(nodes, root, rootcenter, roothalfsize)
    end

    push!(nodes, Node(BoxData(0, Int[], rootcenter, roothalfsize, minlevel), 0, 0, 0))

    firstchildid = _bulkbuildchildren!(
        nodes,
        points,
        buckets,
        rootcenter,
        nsectors,
        childhalfsize,
        minhalfsize,
        0,
        nlevels,
        minvalues,
        maxprotrusion,
        computeprotrusion,
        minlevel + 1,
        root,
        1,
        rootoffset,
    )
    nodes[1] = Node(
        BoxData(0, Int[], rootcenter, roothalfsize, minlevel), 0, 0, firstchildid
    )

    return _finishbulktwontree(nodes, root, rootcenter, roothalfsize)
end

"""
    _finishbulktwontree(nodes, root, rootcenter, roothalfsize)

Materialize a bulk-built node vector as a `TwoNTree`.
"""
function _finishbulktwontree(nodes, root, rootcenter, roothalfsize)
    return TwoNTree(nodes, root, rootcenter, roothalfsize, Vector{Vector{Int}}())
end

# Builds every child of a cell (given its already-computed sector buckets/centers) in Hilbert
# sibling order, wiring up `firstchild`/`nextsibling`/`parent`, and returns the id of the first
# child (0 if none) for the caller to patch into its own node.
"""
    _bulkbuildchildren!(nodes, points, buckets, parentcenter, nsectors, childhalfsize, ...)

Build all occupied children of one cell in Hilbert sibling order.

Returns the global id of the first child, or `0` when no child was built.
"""
function _bulkbuildchildren!(
    nodes,
    points,
    buckets,
    parentcenter,
    nsectors,
    childhalfsize,
    minhalfsize,
    parentdepth,
    nlevels,
    minvalues,
    maxprotrusion,
    computeprotrusion,
    childlevel,
    parentid,
    parentsfcstate,
    rootoffset,
)
    # Occupied sectors in Hilbert sibling order, without heap-allocating a `Vector{Int}` just to
    # sort it once: `nsectors <= 8` always (`N <= 3`), so an `MVector`-backed insertion sort is
    # both stack-allocated and, at this size, no slower than `sort`.
    occupied = MVector{8,Int}(undef)
    positions = MVector{8,Int}(undef)
    n = 0
    hilbertpositions = hilbert_positions[parentsfcstate]
    for idx in 1:nsectors
        isassigned(buckets, idx) || continue
        n += 1
        occupied[n] = idx - 1
        positions[n] = hilbertpositions[idx]
    end
    for i in 2:n
        sec, pos = occupied[i], positions[i]
        j = i - 1
        while j >= 1 && positions[j] > pos
            occupied[j + 1] = occupied[j]
            positions[j + 1] = positions[j]
            j -= 1
        end
        occupied[j + 1] = sec
        positions[j + 1] = pos
    end

    firstchildid = 0
    previd = 0
    for i in 1:n
        sec = occupied[i]
        childsfcstate = hilbert_states[parentsfcstate][sec + 1] + 1
        childid = _bulkbuildnode!(
            nodes,
            points,
            buckets[sec + 1],
            sectorcenter(sec, parentcenter, childhalfsize),
            childhalfsize,
            minhalfsize,
            parentdepth + 1,
            nlevels,
            minvalues,
            maxprotrusion,
            computeprotrusion,
            childlevel,
            parentid,
            sec,
            childsfcstate,
            rootoffset,
        )
        if previd == 0
            firstchildid = childid
        else
            prevvecindex = previd - rootoffset
            prevnode = nodes[prevvecindex]
            nodes[prevvecindex] = Node(
                prevnode.data, childid, prevnode.parent, prevnode.firstchild
            )
        end
        previd = childid
    end
    return firstchildid
end

# `depth`/`level` count from the root; `ids`/`center`/`halfsize` describe this cell; `sector` is
# this cell's sector relative to its own parent (stored on its `BoxData`, needed by
# `locatepoint` if the tree is ever queried after construction).
"""
    _bulkbuildnode!(nodes, points, ids, center, halfsize, ...)

Append one bulk-built node and recursively build its descendants when the cell
does not satisfy the stopping criteria.
"""
function _bulkbuildnode!(
    nodes,
    points,
    ids,
    center,
    halfsize,
    minhalfsize,
    depth,
    nlevels,
    minvalues,
    maxprotrusion,
    computeprotrusion,
    level,
    parentid,
    sector,
    sfcstate,
    rootoffset,
)
    thisvecindex = lastindex(nodes) + 1
    thisid = thisvecindex + rootoffset

    if halfsize <= minhalfsize || depth == nlevels || length(ids) <= minvalues
        push!(nodes, Node(BoxData(sector, ids, center, halfsize, level), 0, parentid, 0))
        return thisid
    end

    nsectors = 1 << length(center)
    childhalfsize = halfsize / 2
    buckets = _bucketbysector(points, ids, center, halfsize, nsectors)

    if !isnan(maxprotrusion)
        for idx in 1:nsectors
            isassigned(buckets, idx) || continue
            if _anyprotrudes(
                points,
                buckets[idx],
                sectorcenter(idx - 1, center, childhalfsize),
                childhalfsize,
                maxprotrusion,
                computeprotrusion,
            )
                push!(
                    nodes,
                    Node(BoxData(sector, ids, center, halfsize, level), 0, parentid, 0),
                )
                return thisid
            end
        end
    end

    push!(nodes, Node(BoxData(sector, Int[], center, halfsize, level), 0, parentid, 0))

    firstchildid = _bulkbuildchildren!(
        nodes,
        points,
        buckets,
        center,
        nsectors,
        childhalfsize,
        minhalfsize,
        depth,
        nlevels,
        minvalues,
        maxprotrusion,
        computeprotrusion,
        level + 1,
        thisid,
        sfcstate,
        rootoffset,
    )
    nodes[thisvecindex] = Node(
        BoxData(sector, Int[], center, halfsize, level), 0, parentid, firstchildid
    )

    return thisid
end

"""
    numberoflevels(halfsize, minhalfsize)

Return the number of bisection levels needed to reach `minhalfsize`.
"""
function numberoflevels(halfsize, minhalfsize)
    return ceil(Int, log2(halfsize / minhalfsize))
end

"""
    roothalfsize(rootsize, minhalfsize)

Return a root halfsize aligned to an integer number of bisections of
`minhalfsize`.
"""
function roothalfsize(rootsize::T, minhalfsize::T) where {T}
    iszero(minhalfsize) && return rootsize
    nlevels = numberoflevels(rootsize, minhalfsize)
    return minhalfsize * T(2)^(nlevels)
end

"""
    cornerpoints(tree::TwoNTree{N,D,T}, node::Int, i)

Return corner `i` of `node`'s box.

Corners are numbered from `1` to `2^N`.
"""
function cornerpoints(tree::TwoNTree{N,D,T}, node::Int, i) where {N,D,T}
    ds = reverse(digits(i - 1; base=2, pad=N))
    for i in eachindex(ds)
        ds[i] == 0 && (ds[i] = -1)
    end

    return center(tree, node) + SVector{N}(ds .* halfsize(tree, node))
end

"""
    treetrait(::Type{<:TwoNTree})

Return [`isTwoNTree`](@ref).
"""
H2Trees.treetrait(::Type{<:TwoNTree}) = isTwoNTree()

"""
    locatepoint(tree, point, targetlevel)

Find the node in `tree` at the given `level` whose box contains `point`.

Starting from the root, descends through child sectors until `level` is reached.
Throws when the existing tree topology has no matching child on the path.
"""
function locatepoint(tree, point, targetlevel)
    tempnode = root(tree)
    for _ in minimumlevel(tree):targetlevel
        if level(tree, tempnode) == targetlevel
            return tempnode
        end

        newsector, _, _ = sectorcentersize(
            point, center(tree, tempnode), halfsize(tree, tempnode)
        )

        for child in children(tree, tempnode)
            if sector(tree, child) == newsector
                tempnode = child
                break
            end
        end
    end

    return error("No suitable node at level $targetlevel found for point $point")
end
