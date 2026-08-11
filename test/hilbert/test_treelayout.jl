using Test
using H2Trees
using StaticArrays
using Random

const HO = H2Trees.HilbertOrdering

# Integration tests: the standalone Hilbert tests prove the curve is correct, these prove the
# final *tree* numbering actually follows it and that renumbering preserved topology/values.

"""
Root-to-node sector path for every node, recovered by walking the tree. The path determines a
node's dyadic grid cell independently of its id, which is what lets the Hilbert order be checked
without trusting the ids under test.
"""
function sectorpaths(tree)
    paths = Dict{Int,Vector{Int}}()
    function visit(node, path)
        paths[node] = copy(path)
        for child in H2Trees.children(tree, node)
            push!(path, H2Trees.sector(tree, child))
            visit(child, path)
            pop!(path)
        end
    end
    visit(H2Trees.root(tree), Int[])
    return paths
end

"""
Integer grid coordinates at the node's own level, from its sector path.
"""
function gridcoordinates(path, ::Val{N}) where {N}
    return ntuple(N) do d
        coordinate = 0
        for sector in path
            coordinate = (coordinate << 1) | ((sector >> (d - 1)) & 1)
        end
        return coordinate
    end
end

function checklayout(tree, ::Val{N}) where {N}
    treelevels = collect(H2Trees.levels(tree))
    rootid = H2Trees.root(tree)
    nnodes = H2Trees.numberofnodes(tree)

    @testset "ids span root:root+n-1 exactly once" begin
        allids = sort(vcat([H2Trees.nodesatlevel(tree, l) for l in treelevels]...))
        @test allids == collect(rootid:(rootid + nnodes - 1))
    end

    @testset "level-major layout" begin
        for level in treelevels
            ids = H2Trees.nodesatlevel(tree, level)
            @test !isempty(ids)
            @test ids == collect(first(ids):last(ids))
        end
        # Adjacent levels touch in id space, so the whole tree is one run of levels.
        for i in 1:(length(treelevels) - 1)
            @test last(H2Trees.nodesatlevel(tree, treelevels[i])) + 1 ==
                first(H2Trees.nodesatlevel(tree, treelevels[i + 1]))
        end
        # The root keeps the requested id and sits first.
        @test H2Trees.nodesatlevel(tree, first(treelevels)) == [rootid]
    end

    paths = sectorpaths(tree)

    @testset "Hilbert order within each level" begin
        for level in treelevels
            depth = level - first(treelevels)
            depth == 0 && continue
            ids = H2Trees.nodesatlevel(tree, level)
            indices = [
                HO.hilbertindex(Val(N), gridcoordinates(paths[id], Val(N)), depth) for
                id in ids
            ]
            # Strictly increasing: increasing node id follows increasing Hilbert index, and no
            # two distinct nodes on a level can share a cell.
            @test issorted(indices)
            @test allunique(indices)
        end
    end

    @testset "topology survived renumbering" begin
        reached = Set{Int}()
        for node in H2Trees.DepthFirstIterator(tree, rootid)
            push!(reached, node)
        end
        @test length(reached) == nnodes

        for level in treelevels, node in H2Trees.nodesatlevel(tree, level)
            childlist = collect(H2Trees.children(tree, node))
            sectors = [H2Trees.sector(tree, c) for c in childlist]
            @test allunique(sectors)

            for child in childlist
                @test H2Trees.parent(tree, child) == node
                @test H2Trees.level(tree, child) == H2Trees.level(tree, node) + 1
                @test H2Trees.halfsize(tree, child) ≈ H2Trees.halfsize(tree, node) / 2
                @test all(
                    H2Trees.center(tree, child) .≈ H2Trees.sectorcenter(
                        H2Trees.sector(tree, child),
                        H2Trees.center(tree, node),
                        H2Trees.halfsize(tree, child),
                    ),
                )
            end

            # firstchild/nextsibling were rewired to Hilbert order, so on a single level (where
            # ids are Hilbert-ordered) the children come out with increasing ids.
            @test issorted(childlist)
            if !isempty(childlist)
                @test H2Trees.firstchild(tree, node) == first(childlist)
                @test H2Trees.nextsibling(tree, last(childlist)) == 0
            end

            node == rootid || @test H2Trees.parent(tree, node) != 0
        end
    end
end

function checkvalues(tree, points, ::Val{N}) where {N}
    @testset "values preserved" begin
        @test sort(H2Trees.values(tree, H2Trees.root(tree))) == collect(eachindex(points))
        # Every point must sit in the leaf its geometry selects.
        maxlevel = maximum(H2Trees.levels(tree))
        for (i, point) in enumerate(points)
            leaf = H2Trees.locatepoint(tree, point, H2Trees.level(tree, findleaf(tree, i)))
            @test i in H2Trees.values(H2Trees.data(tree, leaf))
        end
        @test maxlevel == maximum(H2Trees.levels(tree))
    end
end

findleaf(tree, value) = H2Trees.findleafnode(tree, value)

@testset "tree layout N=$N" for N in 1:3
    Random.seed!(20 + N)

    @testset "uniform random cloud" begin
        points = [SVector(ntuple(_ -> rand(), N)...) for _ in 1:200]
        tree = buildtree(points; builder=TwoNTreeBuilder(; minhalfsize=0.15, minvalues=0))
        checklayout(tree, Val(N))
        checkvalues(tree, points, Val(N))
    end

    @testset "complete grid (balanced, every cell occupied)" begin
        side = 4
        points = SVector{N,Float64}[]
        for idx in Iterators.product(ntuple(_ -> 0:(side - 1), N)...)
            push!(points, SVector(ntuple(d -> idx[d] + 0.5, N)...))
        end
        tree = buildtree(points; builder=TwoNTreeBuilder(; minhalfsize=0.5, minvalues=0))
        checklayout(tree, Val(N))
        checkvalues(tree, points, Val(N))

        @test H2Trees.checkbalancedtree(tree)
        leafids = H2Trees.leaves(tree)
        @test issorted(leafids)
        @test leafids == collect(first(leafids):last(leafids))
        @test leafids == H2Trees.nodesatlevel(tree, maximum(H2Trees.levels(tree)))

        # The downstream chunking guarantee this refactor exists to provide.
        for chunksize in (1, 2, 4, 8, 16, 32)
            chunks = collect(Iterators.partition(first(leafids):last(leafids), chunksize))
            @test vcat(chunks...) == leafids
        end

        # On a complete level, consecutive ids really are face-adjacent boxes (this is exactly
        # the property that fails for sparse occupancy, see below).
        paths = sectorpaths(tree)
        deepest = maximum(H2Trees.levels(tree))
        depth = deepest - minimum(H2Trees.levels(tree))
        coords = [gridcoordinates(paths[id], Val(N)) for id in leafids]
        for i in 1:(length(coords) - 1)
            @test sum(abs.(coords[i + 1] .- coords[i])) == 1
        end
        @test length(leafids) == 1 << (N * depth)
    end

    @testset "sparse balanced grid" begin
        side = 4
        points = SVector{N,Float64}[]
        for idx in Iterators.product(ntuple(_ -> 0:(side - 1), N)...)
            sum(idx) % 3 == 0 || continue
            push!(points, SVector(ntuple(d -> idx[d] + 0.5, N)...))
        end
        tree = buildtree(points; builder=TwoNTreeBuilder(; minhalfsize=0.5, minvalues=0))
        checklayout(tree, Val(N))
        checkvalues(tree, points, Val(N))

        leafids = H2Trees.leaves(tree)
        if H2Trees.checkbalancedtree(tree)
            @test leafids == collect(first(leafids):last(leafids))
        end
        # Deliberately NOT asserted: that consecutive stored leaves are face-adjacent. Only
        # occupied cells are stored, so the stored order is a *subsequence* of the full Hilbert
        # traversal and may jump. `checklayout` pins the surviving guarantee: that
        # the subsequence is in Hilbert order.
        @test allunique(leafids)
    end

    @testset "non-default root id" begin
        points = [SVector(ntuple(_ -> rand(), N)...) for _ in 1:120]
        rootid = 7
        tree = buildtree(
            points; builder=TwoNTreeBuilder(; minhalfsize=0.2, minvalues=0, root=rootid)
        )
        @test H2Trees.root(tree) == rootid
        checklayout(tree, Val(N))
        checkvalues(tree, points, Val(N))
    end

    @testset "single node" begin
        points = [SVector(ntuple(_ -> 0.5, N)...)]
        tree = buildtree(points; builder=TwoNTreeBuilder(; minhalfsize=0.0, minvalues=10))
        @test H2Trees.numberofnodes(tree) == 1
        checklayout(tree, Val(N))
    end
end

@testset "adaptive tree keeps level-major without claiming contiguous leaves" begin
    Random.seed!(3)
    # A dense corner cluster plus a few far points forces leaves onto several levels.
    points = SVector{3,Float64}[]
    for _ in 1:400
        push!(points, SVector(0.02 * rand(), 0.02 * rand(), 0.02 * rand()))
    end
    for _ in 1:8
        push!(points, SVector(0.5 + 0.5 * rand(), 0.5 + 0.5 * rand(), 0.5 + 0.5 * rand()))
    end
    tree = buildtree(points; builder=TwoNTreeBuilder(; minhalfsize=0.0, minvalues=5))

    @test !H2Trees.checkbalancedtree(tree)
    leaflevels = unique(H2Trees.level(tree, leaf) for leaf in H2Trees.leaves(tree))
    @test length(leaflevels) > 1

    # The invariant that DOES survive: every level is still one contiguous Hilbert-ordered block.
    checklayout(tree, Val(3))
    checkvalues(tree, points, Val(3))

    # Deliberately NOT asserted: leaves(tree) == first:last. With leaves on several levels it is
    # impossible to have both level-major contiguity and one contiguous leaf block; level-major
    # is the documented choice.
end

@testset "final ids do not depend on input point order" begin
    Random.seed!(11)
    points = [SVector(rand(), rand(), rand()) for _ in 1:150]
    reference = buildtree(points; builder=TwoNTreeBuilder(; minhalfsize=0.2, minvalues=0))

    for _ in 1:3
        perm = randperm(length(points))
        shuffled = buildtree(
            points[perm]; builder=TwoNTreeBuilder(; minhalfsize=0.2, minvalues=0)
        )

        @test H2Trees.numberofnodes(shuffled) == H2Trees.numberofnodes(reference)
        @test collect(H2Trees.levels(shuffled)) == collect(H2Trees.levels(reference))

        # Geometry keyed by final node id must match exactly: ids are determined by geometry and
        # Hilbert order, not by the order points were processed. (Stored value ids legitimately
        # differ, since the input itself was renumbered, so they are compared through `perm`.)
        for level in H2Trees.levels(reference)
            ids = H2Trees.nodesatlevel(reference, level)
            @test H2Trees.nodesatlevel(shuffled, level) == ids
            for id in ids
                @test H2Trees.sector(shuffled, id) == H2Trees.sector(reference, id)
                @test H2Trees.level(shuffled, id) == H2Trees.level(reference, id)
                @test H2Trees.halfsize(shuffled, id) == H2Trees.halfsize(reference, id)
                @test all(H2Trees.center(shuffled, id) .== H2Trees.center(reference, id))
                @test H2Trees.parent(shuffled, id) == H2Trees.parent(reference, id)
                @test H2Trees.firstchild(shuffled, id) == H2Trees.firstchild(reference, id)
                @test H2Trees.nextsibling(shuffled, id) ==
                    H2Trees.nextsibling(reference, id)
                # Same geometric cell holds the same underlying points, modulo the renumbering.
                @test sort(perm[H2Trees.values(H2Trees.data(shuffled, id))]) ==
                    sort(H2Trees.values(H2Trees.data(reference, id)))
            end
        end
    end
end
