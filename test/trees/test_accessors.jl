module TestAccessors

# Accessor branches that the ordinary tree-building paths never take.
#
# `valuesatnodes` documents a value-to-leaves mapping, but every tree the suite builds stores
# each value in exactly one leaf and numbers values contiguously from one, so the
# "value in several leaves" and "value in no leaf" branches never ran. `leaves(tree, node)`
# likewise has a cached fast path for the root that shadowed its general descent.

using H2Trees
using StaticArrays
using Test

# A hand-built two-level ball tree: root (1) with leaves 2 and 3. Value 2 is stored in BOTH
# leaves, and value 4 in neither, which is what the two uncovered branches are about.
function sharedvaluetree()
    sv(x, y, z) = SVector{3,Float64}(x, y, z)
    nodes = [
        H2Trees.Node(H2Trees.BoundingBallData(Int[], sv(0, 0, 0), 2.0, 1), 0, 0, 2),
        H2Trees.Node(H2Trees.BoundingBallData([1, 2], sv(-1, 0, 0), 1.0, 2), 3, 1, 0),
        H2Trees.Node(H2Trees.BoundingBallData([2, 3], sv(1, 0, 0), 1.0, 2), 0, 1, 0),
    ]
    return H2Trees.BoundingBallTree(nodes, 1, sv(0, 0, 0), 2.0, [[1], [2, 3]])
end

@testset "valuesatnodes maps a value to every leaf holding it" begin
    tree = sharedvaluetree()
    # `numberofvalues` is passed explicitly: the derived count would be 4 (two leaves of two),
    # which happens to match here, but the mapping must be driven by the requested id range,
    # not by how many slots the leaves fill.
    mapping = H2Trees.valuesatnodes(tree; numberofvalues=4)
    @test length(mapping) == 4
    @test mapping[1] == [2]
    # Both leaves, and in ascending node order rather than discovery order.
    @test mapping[2] == [2, 3]
    @test issorted(mapping[2])
    @test mapping[3] == [3]
    # An id no leaf stores must come back as an empty vector, not as an undefined slot.
    @test mapping[4] == Int[]
    @test all(isassigned(mapping, i) for i in eachindex(mapping))
end

@testset "nodesatvalues groups ids by identical leaf sets" begin
    tree = sharedvaluetree()
    grouped = H2Trees.nodesatvalues(tree, H2Trees.valuesatnodes(tree; numberofvalues=4))
    @test grouped[[2]] == [1]
    @test grouped[[2, 3]] == [2]
    @test grouped[[3]] == [3]
    @test grouped[Int[]] == [4]
    # Every value id is accounted for exactly once.
    @test sort(reduce(vcat, values(grouped))) == 1:4
end

@testset "leaves(tree, node) descends from an arbitrary node" begin
    points = [SVector(0.25 * i, 0.25 * j, 0.25 * k) for i in 0:3 for j in 0:3 for k in 0:3]
    tree = H2Trees.buildtree(points; builder=H2Trees.TwoNTreeBuilder(; minvalues=2))
    rootid = H2Trees.root(tree)

    # The root goes through the cached `TreeIndex` list; anything else descends.
    @test H2Trees.leaves(tree) == H2Trees.leaves(tree, rootid)
    # A copy, so a caller mutating the result cannot corrupt the cached index.
    cached = H2Trees.leaves(tree)
    push!(cached, -1)
    @test H2Trees.leaves(tree) != cached

    internal = first(
        Iterators.filter(
            n -> !H2Trees.isleaf(tree, n) && n != rootid, H2Trees.depthfirstnodes(tree)
        ),
    )
    sub = H2Trees.leaves(tree, internal)
    @test !isempty(sub)
    @test all(H2Trees.isleaf(tree, n) for n in sub)
    # Every leaf below `internal` and no others: check against an independent upward walk.
    expected = filter(H2Trees.leaves(tree)) do leaf
        internal in H2Trees.ParentUpwardsIterator(tree, leaf)
    end
    @test sort(sub) == sort(expected)

    # A leaf is its own only leaf.
    someleaf = first(H2Trees.leaves(tree))
    @test H2Trees.leaves(tree, someleaf) == [someleaf]
end

@testset "istranslatingnode is false where nothing is well separated" begin
    points = [SVector(0.25 * i, 0.25 * j, 0.25 * k) for i in 0:3 for j in 0:3 for k in 0:3]
    tree = H2Trees.buildtree(points; builder=H2Trees.TwoNTreeBuilder(; minvalues=2))

    # The root has no parent, so nothing can be well separated from it.
    @test !H2Trees.istranslatingnode(tree, H2Trees.root(tree))
    # `0` is the "no such node" sentinel and must answer without touching the tree.
    @test !H2Trees.istranslatingnode(tree, 0)

    # And it must agree with the iterator it is a shortcut for, on every node.
    for node in H2Trees.depthfirstnodes(tree)
        @test H2Trees.istranslatingnode(tree, node) ==
            !isempty(collect(H2Trees.TranslatingNodesIterator(tree, node)))
    end
end

@testset "mintranslationlevel falls back to the finest level" begin
    # A tree small enough that no pair is ever well separated: the documented fallback is the
    # last level rather than an error or a wrong `0`.
    points = [SVector(0.01 * i, 0.0, 0.0) for i in 0:9]
    tree = H2Trees.buildtree(points; builder=H2Trees.TwoNTreeBuilder(; minvalues=100))
    @test all(!H2Trees.istranslatingnode(tree, n) for n in H2Trees.depthfirstnodes(tree))
    @test H2Trees.mintranslationlevel(tree) == H2Trees.levels(tree)[end]
end

end # module TestAccessors
