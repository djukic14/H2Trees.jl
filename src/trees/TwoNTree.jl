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
    _uniformseparationdepth(points, indices, center, halfsize, cutoff::Int=1)

Number of bisection levels needed until each cell holds at most `cutoff` points.

This uses the same sector rule as `TwoNTree` without building a tree. `cutoff` is the caller's
`minvalues`, floored at 1.

# Coincident points

Coincident points cannot be separated by bisection: they stay in the same sector at every level.

This throws `ArgumentError` when that would prevent termination, i.e. when a cell still holds more
than `cutoff` points after `halfsize` has reached `0.0`. It is NOT a general duplicate check:
with `cutoff > 1`, coincident points that fit inside one `<= cutoff` cell simply stop the
recursion normally and are never noticed. Duplicate input points are not permitted and are not
validated; see [`bulkbuildtwontree`](@ref).
"""
function _uniformseparationdepth(points, indices, center, halfsize, cutoff::Int=1)
    npoints = length(indices)
    npoints <= max(cutoff, 1) && return 0

    # Scratch is allocated once and partitioned in place by sector.
    workingindices = collect(Int, indices)
    return _separationdepth!(
        points,
        workingindices,
        similar(workingindices),
        Vector{Int}(undef, npoints),
        1,
        npoints,
        center,
        halfsize,
        max(cutoff, 1),
    )
end

# Recursive worker for `_uniformseparationdepth`. `indices[lo:hi]` is this cell's point ids;
# `scratch`/`sectors` are shared buffers reused at every depth. Returns the number of further
# bisections this cell needs.
function _separationdepth!(
    points, indices, scratch, sectors, lo::Int, hi::Int, center, halfsize, cutoff::Int
)
    (hi - lo + 1) <= cutoff && return 0
    if iszero(halfsize)
        return throw(
            ArgumentError(
                "cannot separate $(hi - lo + 1) coincident points (or points closer than " *
                "floating-point precision can resolve) into distinct boxes; deduplicate the " *
                "input points before building a tree",
            ),
        )
    end

    # `sc` (from `sectorcentersize`) always lands in `0:2^N-1`, and `N <= 3` for every supported
    # dimension, so per-sector counters fit in fixed stack storage.
    nsectors = 1 << length(center)
    counts = zero(MVector{8,Int})
    @inbounds for k in lo:hi
        sector, _, _ = sectorcentersize(points[indices[k]], center, halfsize)
        sectors[k] = sector
        counts[sector + 1] += 1
    end

    # Counting-sort offsets: sector `s` occupies `offsets[s] : offsets[s] + counts[s] - 1`.
    offsets = zero(MVector{8,Int})
    next = lo
    @inbounds for sector in 1:nsectors
        offsets[sector] = next
        next += counts[sector]
    end

    cursor = copy(offsets)
    @inbounds for k in lo:hi
        slot = sectors[k] + 1
        scratch[cursor[slot]] = indices[k]
        cursor[slot] += 1
    end
    @inbounds for k in lo:hi
        indices[k] = scratch[k]
    end

    # This cell performs a real split, so it counts as one level even when the split already
    # separates everything (the recursive calls then contribute 0).
    childhalfsize = halfsize / 2
    maxchilddepth = 0
    @inbounds for sector in 1:nsectors
        iszero(counts[sector]) && continue
        maxchilddepth = max(
            maxchilddepth,
            _separationdepth!(
                points,
                indices,
                scratch,
                sectors,
                offsets[sector],
                offsets[sector] + counts[sector] - 1,
                sectorcenter(sector - 1, center, childhalfsize),
                childhalfsize,
                cutoff,
            ),
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

  - the general adaptive rule: `length(ids) <= minvalues`, protrusion, `halfsize <= minhalfsize`,
    and the bisection depth [`_uniformseparationdepth`](@ref) computes;
  - `minvalues == 0` with `NoProtrusionCheck`: a uniform-depth rule bounded ONLY by
    `halfsize <= minhalfsize`, with no conformance computation at all (points landing in the same
    cell all the way down to `minhalfsize` simply share a leaf), implemented by passing
    `_bulkbuildnode!`/`_bulkbuildchildren!` a sentinel `nlevels = typemax(Int)` so the
    depth-based stopping condition never fires, skipping `_uniformseparationdepth` entirely
    (see the `nlevels` computation below).

!!! warning "Duplicate points are not permitted and are not validated"

    Two points at the identical location can never be placed in different boxes, however deep
    the tree goes. Passing them is a caller error, and H2Trees does NOT reliably detect it.

    [`_uniformseparationdepth`](@ref) throws an `ArgumentError` only when duplicates would
    actually prevent termination: more than `minvalues` points sharing one location, so no
    amount of subdivision brings a cell below the `minvalues` stopping rule. Fewer duplicates
    than that stop the scan normally and pass unnoticed; they end up sharing a leaf, exactly as
    they already do on the `minvalues == 0` uniform path. Deduplicate before building.

    Note that "identical" here also covers points too close for the coordinate type to resolve
    at the relevant scale, which no exact-equality check would catch either.

Either way the recursion itself is purely geometric: children are appended in ascending sector
order, and the resulting node ids and sibling links are temporary. [`_finishbulktwontree`](@ref)
then assigns the final layout: level-major, Hilbert-ordered within each level, so
`firstchild`/`nextsibling` and node ids alike are deterministic regardless of point processing
order. See [`_hilbertlevelmajornodes`](@ref) for the guarantees that layout provides.

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
    # Checked HERE rather than left to `HilbertOrdering`, which also rejects N outside 1:3 but
    # is only reached after construction: `_separationdepth!` counts sectors into a fixed
    # `MVector{8}` under `@inbounds`, so `N >= 4` (16+ sectors) writes past it and segfaults the
    # process long before any Julia-level error can be raised.
    1 <= N <= 3 || throw(
        ArgumentError(
            "TwoNTree supports 1, 2 and 3 dimensions, got $N. The sector bookkeeping and the " *
            "Hilbert ordering are both built around at most 2^3 children per node.",
        ),
    )
    # `minhalfsize` is the ONLY thing bounding recursion in the `minvalues=0`+`NoProtrusionCheck`
    # uniform-bulk regime below (it deliberately skips `_uniformseparationdepth`'s depth cap);
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
    # separated into distinct leaves; a sentinel nlevels the depth check in
    # `_bulkbuildnode!`/`_bulkbuildchildren!` never reaches, so `_uniformseparationdepth` is never
    # called at all for this regime. See the docstring above.
    nlevels = if iszero(minvalues) && protrusion isa NoProtrusionCheck
        typemax(Int)
    else
        # `minvalues` as the scan's cutoff: the adaptive path stops a branch the moment it holds
        # `<= minvalues` points, so the depth at which every branch reaches that IS the cap. The
        # scan used to recurse to one point per cell, measuring depth the builder can never use.
        # For `minvalues=0` (only reachable here with a real protrusion check) the cutoff floors
        # to 1, i.e. exactly the previous behaviour.
        max(
            1,
            _uniformseparationdepth(
                points, eachindex(points), rootcenter, roothalfsize, minvalues
            ),
        )
    end

    nodes = Vector{Node{BoxData{N,T}}}()
    # Conservative upper bound on final node count: a full `2^N`-ary tree over `length(points)`
    # leaves has at most `(2^N * length(points) - 1) / (2^N - 1)` nodes; the simpler
    # `(2^N) * length(points)` bound used here is looser but cheap to compute and avoids
    # `sizehint!` itself becoming a source of repeated reallocation on the (much more common)
    # shallower trees.
    sizehint!(nodes, min(2 * length(points), 1 + (1 << N) * length(points)))
    rootoffset = root - 1

    # The root is subject to the SAME stop rules as every other node: `_bulkbuildnode!` stops
    # on `halfsize <= minhalfsize || depth == nlevels || length(ids) <= minvalues`, and the root
    # is just that test at depth 0. Checking only `minhalfsize` here meant a tree small enough to
    # stop immediately was still split once, because `minvalues` was not consulted until after
    # the root's children had been created.
    if roothalfsize <= minhalfsize || length(points) <= minvalues
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
    # root; see the docstring above.
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
        rootoffset,
    )
    nodes[1] = Node(
        BoxData(0, Int[], rootcenter, roothalfsize, minlevel), 0, 0, firstchildid
    )

    return _finishbulktwontree(nodes, root, rootcenter, roothalfsize)
end

"""
    _finishbulktwontree(nodes, root, rootcenter, roothalfsize)

Materialize a bulk-built node vector as a `TwoNTree`, assigning final node ids first.

The bulk builder's ids and sibling links are temporary (ascending sector order, depth-first
numbering). This is the one point where they become final: [`_hilbertlevelmajornodes`](@ref)
rewrites both before the `TwoNTree` and [`TreeIndex`](@ref) exist. Nothing outside construction
can observe the temporary ids.

Reordering lives here rather than in a `TwoNTree(nodes, ...)` constructor on purpose: a general
constructor that silently permuted already-materialized nodes would invalidate node ids held by
plans and other callers.
"""
function _finishbulktwontree(
    nodes::Vector{Node{BoxData{N,T}}}, root, rootcenter, roothalfsize
) where {N,T}
    return TwoNTree(
        _hilbertlevelmajornodes(nodes, root, Val(N)),
        root,
        rootcenter,
        roothalfsize,
        Vector{Vector{Int}}(),
    )
end

"""
    _hilbertlevelmajornodes(nodes, root, ::Val{N})

Return a new node vector renumbered level-major with Hilbert order inside each level.

Walks the temporary topology depth-first with children visited in Hilbert order, propagating the
orientation state from the root via `HilbertOrdering`. A preorder walk in Hilbert child order
visits each level's nodes in exactly that level's Hilbert order: a node's Hilbert index is
`(parent index) * 2^N + (position within parent)`, so grouping by ancestor in ancestor order and
ordering within a group by position *is* Hilbert order, so bucketing nodes by level during the
walk needs no sorting or per-node Hilbert keys.

Sibling links are rebuilt too, not just permuted: `firstchild`/`nextsibling` define the child
traversal order that `children` exposes, so they must follow the Hilbert order the ids do.

The walk is iterative rather than recursive so this adds no stack-depth limit of its own beyond
the one the recursive builder already imposes. Runs in `O(length(nodes))`; all scratch storage is
released when construction finishes, and no Hilbert state is stored on the returned nodes.

The returned vector is allocated at exactly `length(nodes)`, which also drops the oversized
backing buffer the bulk builder's `sizehint!` leaves behind (that hint is a loose upper bound on
the node count, so a finished tree could otherwise retain several megabytes of unused capacity
for its whole lifetime.
"""
function _hilbertlevelmajornodes(
    nodes::Vector{Node{BoxData{N,T}}}, root::Int, ::Val{N}
) where {N,T}
    nnodes = length(nodes)
    # A lone root is already trivially level-major and Hilbert-ordered, but still gets copied
    # into an exactly-sized vector for the same reason the general path below builds one: the
    # bulk builder `sizehint!`s `nodes` to a conservative upper bound on the node count, and
    # that oversized backing buffer would otherwise be retained by the tree for its whole
    # lifetime (worst case here: a 1-node tree holding a buffer sized for `2 * npoints`).
    nnodes <= 1 && return copy(nodes)

    rootoffset = root - 1
    nsectors = 1 << N

    minlevel = level(first(nodes).data)
    maxlevel = minlevel
    for node in nodes
        nodelevel = level(node.data)
        nodelevel > maxlevel && (maxlevel = nodelevel)
    end

    nodesbylevel = [Int[] for _ in minlevel:maxlevel]
    # Rebuilt sibling topology, indexed by OLD vector index, filled during the walk below.
    hilbertfirstchild = zeros(Int, nnodes)
    hilbertnextsibling = zeros(Int, nnodes)

    # `nsectors <= 8` for every supported dimension, so children of one node fit in fixed
    # stack storage and need no per-node heap allocation.
    childids = MVector{8,Int}(undef)
    childpositions = MVector{8,Int}(undef)

    stack = [(root, HilbertOrdering.initialstate(Val(N)))]
    while !isempty(stack)
        nodeid, state = pop!(stack)
        vecindex = nodeid - rootoffset
        node = nodes[vecindex]
        push!(nodesbylevel[level(node.data) - minlevel + 1], nodeid)

        nchildren = 0
        childid = node.firstchild
        while !iszero(childid)
            childnode = nodes[childid - rootoffset]
            nchildren += 1
            childids[nchildren] = childid
            childpositions[nchildren] = HilbertOrdering.hilbertposition(
                Val(N), state, sector(childnode.data)
            )
            childid = childnode.nextsibling
        end

        # Insertion sort: at most 8 elements, so this beats a general sort and stays allocation
        # free on the fixed-size scratch above.
        for i in 2:nchildren
            id, position = childids[i], childpositions[i]
            j = i - 1
            while j >= 1 && childpositions[j] > position
                childids[j + 1] = childids[j]
                childpositions[j + 1] = childpositions[j]
                j -= 1
            end
            childids[j + 1] = id
            childpositions[j + 1] = position
        end

        if nchildren > 0
            hilbertfirstchild[vecindex] = childids[1]
            for i in 1:(nchildren - 1)
                hilbertnextsibling[childids[i] - rootoffset] = childids[i + 1]
            end
            hilbertnextsibling[childids[nchildren] - rootoffset] = 0
        end

        # Pushed in reverse so the stack pops them in Hilbert order.
        for i in nchildren:-1:1
            id = childids[i]
            push!(
                stack,
                (
                    id,
                    HilbertOrdering.hilbertnextstate(
                        Val(N), state, sector(nodes[id - rootoffset].data)
                    ),
                ),
            )
        end
    end

    newtoold = Vector{Int}(undef, nnodes)
    next = 1
    for levelnodes in nodesbylevel
        for oldid in levelnodes
            newtoold[next] = oldid
            next += 1
        end
    end
    # A walk that did not reach every node would mean the temporary topology was malformed;
    # catching it here beats a confusing `undef` id surfacing later from `TreeIndex`.
    next == nnodes + 1 || error(
        "Hilbert reorder reached $(next - 1) of $nnodes nodes; bulk-built topology is not a " *
        "single tree rooted at $root",
    )

    oldtonew = Vector{Int}(undef, nnodes)
    for (newindex, oldid) in enumerate(newtoold)
        oldtonew[oldid - rootoffset] = newindex + rootoffset
    end
    # `0` is the "no such node" sentinel used by parent/firstchild/nextsibling and must survive
    # remapping as itself; global ids are otherwise offset by `root`, never equal to the vector
    # index unless `root == 1`.
    remap(id::Int) = iszero(id) ? 0 : oldtonew[id - rootoffset]

    newnodes = Vector{Node{BoxData{N,T}}}(undef, nnodes)
    for (newindex, oldid) in enumerate(newtoold)
        oldvecindex = oldid - rootoffset
        oldnode = nodes[oldvecindex]
        # Only identity and topology references change; `data` (sector, values, center,
        # halfsize, level) is carried over untouched.
        newnodes[newindex] = Node(
            oldnode.data,
            remap(hilbertnextsibling[oldvecindex]),
            remap(oldnode.parent),
            remap(hilbertfirstchild[oldvecindex]),
        )
    end

    return newnodes
end

# Builds every child of a cell (given its already-computed sector buckets/centers) in plain
# ascending sector order, wiring up `firstchild`/`nextsibling`/`parent`, and returns the id of
# the first child (0 if none) for the caller to patch into its own node.
"""
    _bulkbuildchildren!(nodes, points, buckets, parentcenter, nsectors, childhalfsize, ...)

Build all occupied children of one cell in ascending sector order.

Returns the global id of the first child, or `0` when no child was built.

Sibling order here is deliberately *temporary*. Node ids and sibling links produced by the
recursive bulk build are construction details; [`_finishbulktwontree`](@ref) rewrites both into
the final level-major, Hilbert-within-level layout before any `TwoNTree`/`TreeIndex` exists.
Keeping Hilbert state out of this recursion leaves it a pure geometric partitioner.
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
    rootoffset,
)
    firstchildid = 0
    previd = 0
    for idx in 1:nsectors
        isassigned(buckets, idx) || continue
        sec = idx - 1
        childid = _bulkbuildnode!(
            nodes,
            points,
            buckets[idx],
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
