using Test
using CompScienceMeshes
using StaticArrays
using H2Trees
using BEAST

@testset "Empty tree" begin
    tree = TwoNTree(SVector(0.0, 0.0, 0.0), 1.0)

    @test H2Trees.root(tree) == 1
    @test H2Trees.center(tree, 1) == SVector(0.0, 0.0, 0.0)
    @test H2Trees.halfsize(tree, 1) == 1.0
    @test H2Trees.level(tree, 1) == 1
    @test H2Trees.isin(tree, 1, SVector(1.0, 1.0, 1.0))
    @test !H2Trees.isin(tree, 1, SVector(1.0 + 1.0e-8, 0.0, 0.0))
    @test !H2Trees.isin(tree, 1, SVector(1.1, 1.1, 1.1))
    @test H2Trees.values(tree, 1) == Int[]
    @test H2Trees.values(tree, H2Trees.leaves(tree)) == Int[]
    @test H2Trees.LevelIterator(tree, 1) == [1]
    @test H2Trees.treeindex(tree).nodes_by_level == [[1]]
    @test H2Trees.treetrait(tree) == H2Trees.isTwoNTree()
end

@testset "Filled Tree" begin
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere4.in")
    )

    points = vertices(m)

    root = 2
    minlevel = 2

    tree = buildtree(
        points;
        builder=TwoNTreeBuilder(;
            minhalfsize=0.1, root=root, minlevel=minlevel, minvalues=10
        ),
    )

    @test H2Trees.halfsizes(tree) == [0.8, 0.4, 0.2, 0.1]

    @test H2Trees.values(tree, H2Trees.root(tree)) ==
        H2Trees.values(tree, H2Trees.leaves(tree, H2Trees.root(tree)))

    appended = [-1]
    @test H2Trees.appendvalues!(appended, tree, H2Trees.root(tree)) === appended
    @test appended == [-1; H2Trees.values(tree, H2Trees.root(tree))]

    nodesappended = [-2]
    @test H2Trees.appendvalues!(
        nodesappended, tree, H2Trees.leaves(tree, H2Trees.root(tree))
    ) === nodesappended
    @test nodesappended == [-2; H2Trees.values(tree, H2Trees.root(tree))]

    visited = Int[]
    @test H2Trees.foreachvalue(tree, H2Trees.root(tree)) do value
        push!(visited, value)
    end === nothing
    @test visited == H2Trees.values(tree, H2Trees.root(tree))
    @test H2Trees.anyvalue(==(first(visited)), tree, H2Trees.root(tree))
    @test !H2Trees.anyvalue(==(0), tree, H2Trees.root(tree))

    for node in H2Trees.DepthFirstIterator(tree)
        @test H2Trees.samelevelnodes(tree, node) ==
            H2Trees.nodesatlevel(tree, H2Trees.level(tree, node))
    end

    for leaf in H2Trees.leaves(tree)
        for point in H2Trees.values(tree, leaf)
            @test H2Trees.isin(tree, leaf, vertices(m)[point])
        end
    end

    valuesatnodes = H2Trees.valuesatnodes(tree)
    @test length(valuesatnodes) == length(points)
    for (functionid, value) in enumerate(valuesatnodes)
        @test length(value) == 1
        @test functionid in H2Trees.values(tree, value[1])
    end
    nodesatvalues = H2Trees.nodesatvalues(tree)
    for (key, value) in nodesatvalues
        @test length(key) == 1
        key = key[1]
        @test sort(value) == sort(H2Trees.values(tree, key))
    end

    maximumlevel = H2Trees.levels(tree)[end] - minlevel + 1

    leaflevels = sort(unique(H2Trees.level.(Ref(tree), H2Trees.leaves(tree))))[2:end]

    nodes = Int[]
    for i in H2Trees.DepthFirstIterator(tree, root)
        @test any(i .== H2Trees.LevelIterator(tree, H2Trees.level(tree, i)))

        for node in H2Trees.SameLevelIterator(tree, i)
            @test H2Trees.level(tree, node) == H2Trees.level(tree, i)
        end

        for node in H2Trees.TranslatingNodesIterator(tree, i)
            @test H2Trees.level(tree, node) == H2Trees.level(tree, i)
            @test !H2Trees.isnear(tree, node, i)
        end

        valuesonlevel = Int[]

        for node in H2Trees.TranslatingNodesIterator(tree, i)
            append!(valuesonlevel, H2Trees.values(tree, node))
        end

        for node in H2Trees.NotTranslatingNodesIterator(tree, i)
            append!(valuesonlevel, H2Trees.values(tree, node))
        end

        if !(H2Trees.level(tree, i) in leaflevels)
            @test sort!(valuesonlevel) == Array(1:length(points))
        end
        push!(nodes, i)

        sector = H2Trees.sector(tree, i)
        oppositesector = H2Trees.oppositesector(tree, i)

        if sector == 0
            @test oppositesector == 7
        elseif sector == 1
            @test oppositesector == 6
        elseif sector == 2
            @test oppositesector == 5
        elseif sector == 3
            @test oppositesector == 4
        elseif sector == 4
            @test oppositesector == 3
        elseif sector == 5
            @test oppositesector == 2
        elseif sector == 6
            @test oppositesector == 1
        elseif sector == 7
            @test oppositesector == 0
        end

        @test H2Trees.levelindex(tree, i) == H2Trees.level(tree, i) - minlevel + 1
        i == root && continue

        pminuschild = H2Trees.parentcenterminuschildcenter(tree, i)

        @test pminuschild ≈
            H2Trees.center(tree, H2Trees.parent(tree, i)) - H2Trees.center(tree, i)

        for cornerid in 1:8
            corner = H2Trees.cornerpoints(tree, i, cornerid)

            if cornerid == 1
                @test corner ≈ H2Trees.center(tree, i) .- H2Trees.halfsize(tree, i)

            elseif cornerid == 2
                @test corner ≈
                    H2Trees.center(tree, i) +
                      H2Trees.halfsize(tree, i) .* SVector(-1, -1, 1)

            elseif cornerid == 3
                @test corner ≈
                    H2Trees.center(tree, i) +
                      H2Trees.halfsize(tree, i) .* SVector(-1, 1, -1)

            elseif cornerid == 4
                @test corner ≈
                    H2Trees.center(tree, i) +
                      H2Trees.halfsize(tree, i) .* SVector(-1, 1, 1)

            elseif cornerid == 5
                @test corner ≈
                    H2Trees.center(tree, i) +
                      H2Trees.halfsize(tree, i) .* SVector(1, -1, -1)
            elseif cornerid == 6
                @test corner ≈
                    H2Trees.center(tree, i) +
                      H2Trees.halfsize(tree, i) .* SVector(1, -1, 1)

            elseif cornerid == 7
                @test corner ≈
                    H2Trees.center(tree, i) +
                      H2Trees.halfsize(tree, i) .* SVector(1, 1, -1)
            else
                cornerid == 8
                @test corner ≈
                    H2Trees.center(tree, i) +
                      H2Trees.halfsize(tree, i) .* SVector(1, 1, 1)
            end
        end
    end

    @test sort!(nodes) == Array(root:(length(tree.nodes) + root - 1))

    for level in H2Trees.levels(tree)
        level in leaflevels && continue

        iteratorvalues = H2Trees.values(tree, H2Trees.LevelIterator(tree, level))

        values = Int[]
        for node in H2Trees.LevelIterator(tree, level)
            append!(values, H2Trees.values(tree, node))
        end
        @test iteratorvalues == values
        @test sort(iteratorvalues) == Array(1:length(points))
        @test sort(values) == Array(1:length(points))
    end

    leaves = H2Trees.leaves(tree)

    for leaf in leaves
        @test H2Trees.isleaf(tree, leaf)
        @test tree(leaf).firstchild == 0
    end

    tree2 = TwoNTree(SVector(0.0, 0.0, 0.0), 0.1; root=root, minlevel=minlevel)
    @test H2Trees.root(tree2) == root
    @test H2Trees.center(tree2, root) == SVector(0.0, 0.0, 0.0)
    @test H2Trees.halfsize(tree2, root) == 0.1
    @test H2Trees.level(tree2, root) == minlevel

    @test H2Trees.treetrait(tree) == H2Trees.isTwoNTree()
    @test H2Trees.treetrait(tree2) == H2Trees.isTwoNTree()

    @test H2Trees.minhalfsize(tree) ≈ H2Trees.halfsize(tree) * (1 / 2)^(maximumlevel - 1)

    for value in eachindex(points)
        @test H2Trees.findleafnode(tree, value) ∈ leaves
        @test value ∈ H2Trees.values(tree, H2Trees.findleafnode(tree, value))
    end

    @test H2Trees.findleafnode(tree, length(points) + 1) == 0
    @test H2Trees.findleafnode(tree, -1) == 0

    leafclusters = H2Trees.leafclusters(tree)
    leavesvalues = H2Trees.values(tree, H2Trees.leaves(tree))
    for (i, leaf) in enumerate(H2Trees.leaves(tree))
        @test leafclusters[i] == H2Trees.values(tree, leaf)
        @test issubset(leafclusters[i], leavesvalues)
    end
end

@testset "Tree builders and index" begin
    points = [
        SVector(0.0, 0.0, 0.0),
        SVector(1.0, 1.0, 1.0),
        SVector(0.25, 0.25, 0.25),
        SVector(0.75, 0.75, 0.75),
    ]

    tree = TwoNTree(
        points;
        builder=TwoNTreeBuilder(;
            minhalfsize=0.25,
            minlevel=AutoMinLevel(),
            minvalues=1,
            protrusion=NoProtrusionCheck(),
        ),
    )

    @test H2Trees.minimumlevel(tree) == 1
    @test H2Trees.treeindex(tree).nodes_by_level == H2Trees.nodesatlevel(tree)
    @test H2Trees.depthfirstnodes(tree) == collect(Int, H2Trees.DepthFirstIterator(tree))
    @test H2Trees.treeindex(tree).leaves == H2Trees.leaves(tree)
    # `leaves(tree)` reuses the cached `TreeIndex.leaves` for the root case, but must hand back an
    # independently owned copy: mutating the result must not corrupt the cache a later call reads.
    @test H2Trees.leaves(tree) !== H2Trees.treeindex(tree).leaves
    let returned = H2Trees.leaves(tree)
        push!(returned, -1)
        @test -1 ∉ H2Trees.treeindex(tree).leaves
        @test -1 ∉ H2Trees.leaves(tree)
    end
    @test H2Trees.levelprotrusions(tree) == H2Trees.maxprotrusion(tree)
    @test H2Trees.protrusionreport(tree).level in H2Trees.levels(tree)
    @test (:nodes, :root, :center, :halfsize, :index) ⊆ propertynames(tree)
    @test :name ∉ propertynames(tree)
    @test :super ∉ propertynames(tree)

    rawtwontreenode = H2Trees.Node(
        H2Trees.BoxData(0, Int[], SVector(0.0, 0.0, 0.0), 1.0, 1), 0, 0, 0
    )
    rawtwontree = TwoNTree([rawtwontreenode], 1, SVector(0.0, 0.0, 0.0), 1.0, [[1]])
    @test H2Trees.treeindex(rawtwontree) isa H2Trees.TreeIndex
    @test H2Trees.nodesatlevel(rawtwontree) == [[1]]
    typedrawtwontree = TwoNTree{3,H2Trees.BoxData{3,Float64},Float64}(
        [rawtwontreenode], 1, SVector(0.0, 0.0, 0.0), 1.0, [[1]]
    )
    @test H2Trees.treeindex(typedrawtwontree) isa H2Trees.TreeIndex
    @test H2Trees.nodesatlevel(typedrawtwontree) == [[1]]

    kmeansbuilder = KMeansTreeBuilder(; numberofclusters=2)
    metisbuilder = MetisTreeBuilder(; numdivisions=2)
    hybridbuilder = SimpleHybridTreeBuilder(; hybridhalfsize=0.25)
    @test kmeansbuilder.numberofclusters == 2
    @test haskey(kmeansbuilder.splitterkwargs, :rng)
    @test metisbuilder.numdivisions == 2
    @test !metisbuilder.splitunconnectedpartitions
    @test hybridbuilder.hybridhalfsize == 0.25
    @test fieldtype(typeof(kmeansbuilder), :minlevel) !== Any
    @test fieldtype(typeof(kmeansbuilder), :balldata) !== Any
    @test fieldtype(typeof(kmeansbuilder), :updateradii) !== Any
    @test fieldtype(typeof(metisbuilder), :minlevel) !== Any
    @test fieldtype(typeof(hybridbuilder), :hybridhalfsize) !== Any

    block = buildtree(
        points,
        points;
        builder=BlockTreeBuilder(;
            test=TwoNTreeBuilder(; minhalfsize=0.25, minvalues=1),
            trial=TwoNTreeBuilder(; minhalfsize=0.25, minvalues=2),
        ),
    )

    @test H2Trees.minhalfsize(H2Trees.testtree(block)) ==
        H2Trees.minhalfsize(H2Trees.trialtree(block))
    @test_throws ArgumentError BlockTreeBuilder(;
        test=TwoNTreeBuilder(; minhalfsize=0.25), trial=TwoNTreeBuilder(; minhalfsize=0.5)
    )
end

@testset "BoundingBallTree splitwrapper compatibility" begin
    points = [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(0.0, 1.0), SVector(1.0, 1.0)]

    oldsplit(points, globalpointids, numsplits) = (
        [globalpointids[1:2], globalpointids[3:4]],
        [SVector(0.5, 0.0), SVector(0.5, 1.0)],
        [0.5, 0.5],
    )
    newsplit(points, globalpointids, level, numsplits) = (
        [globalpointids[1:2], globalpointids[3:4]],
        [SVector(level, 0.0), SVector(level, 1.0)],
        [0.5, 0.5],
    )

    roottree = H2Trees.BoundingBallTree(SVector(0.0, 0.0), 1.0)
    @test H2Trees.nodesatlevel(roottree, H2Trees.minimumlevel(roottree)) == [1]
    @test H2Trees.treeindex(roottree).leaves == [1]
    @test (:nodes, :root, :center, :radius, :index) ⊆ propertynames(roottree)
    @test :name ∉ propertynames(roottree)
    @test :super ∉ propertynames(roottree)

    rawballnode = H2Trees.Node(
        H2Trees.BoundingBallData([1], SVector(0.0, 0.0), 1.0, 1), 0, 0, 0
    )
    rawballtree = H2Trees.BoundingBallTree([rawballnode], 1, SVector(0.0, 0.0), 1.0, [[1]])
    @test H2Trees.treeindex(rawballtree) isa H2Trees.TreeIndex
    @test H2Trees.nodesatlevel(rawballtree) == [[1]]
    typedrawballtree = H2Trees.BoundingBallTree{
        2,H2Trees.BoundingBallData{2,Float64},Float64
    }(
        [rawballnode], 1, SVector(0.0, 0.0), 1.0, [[1]]
    )
    @test H2Trees.treeindex(typedrawballtree) isa H2Trees.TreeIndex
    @test H2Trees.nodesatlevel(typedrawballtree) == [[1]]

    oldtree = H2Trees.buildtree(
        points;
        builder=BoundingBallTreeBuilder(; splitter=oldsplit, numsplits=2, minvalues=1),
    )
    newtree = H2Trees.buildtree(
        points;
        builder=BoundingBallTreeBuilder(; splitter=newsplit, numsplits=2, minvalues=1),
    )
    buildertree = H2Trees.BoundingBallTree(
        points;
        builder=BoundingBallTreeBuilder(;
            splitter=newsplit, numsplits=2, minvalues=1, splitterkwargs=(;)
        ),
    )

    @test H2Trees.nodesatlevel(oldtree, H2Trees.minimumlevel(oldtree)) == [1]
    @test H2Trees.nodesatlevel(newtree, H2Trees.minimumlevel(newtree)) == [1]
    @test H2Trees.nodesatlevel(buildertree, H2Trees.minimumlevel(buildertree)) == [1]
    @test H2Trees.center(newtree, first(H2Trees.nodesatlevel(newtree, 2))) ==
        SVector(1.0, 0.0)
end

@testset "BoundingBallTree balanceleaves!" begin
    nodes = [
        H2Trees.Node(H2Trees.BoundingBallData(Int[], SVector(0.0, 0.0), 2.0, 1), 0, 0, 2),
        H2Trees.Node(H2Trees.BoundingBallData([1], SVector(-1.0, 0.0), 1.0, 2), 3, 1, 0),
        H2Trees.Node(H2Trees.BoundingBallData(Int[], SVector(1.0, 0.0), 1.0, 2), 0, 1, 4),
        H2Trees.Node(H2Trees.BoundingBallData([2], SVector(1.0, 0.0), 0.5, 3), 0, 3, 0),
    ]
    tree = H2Trees.BoundingBallTree(nodes, 1, SVector(0.0, 0.0), 2.0, [[1], [2, 3], [4]])

    @test sort(H2Trees.leaves(tree)) == [2, 4]
    @test H2Trees.values(tree, H2Trees.root(tree)) == [1, 2]
    @test H2Trees.balanceleaves!(tree) === tree
    @test sort(H2Trees.leaves(tree)) == [4, 5]
    @test all(H2Trees.level(tree, leaf) == 3 for leaf in H2Trees.leaves(tree))
    @test H2Trees.firstchild(tree, 2) == 5
    @test H2Trees.parent(tree, 5) == 2
    @test H2Trees.data(tree, 2).values == Int[]
    @test H2Trees.values(tree, 5) == [1]
    @test H2Trees.values(tree, 4) == [2]
    @test sort(H2Trees.values(tree, H2Trees.root(tree))) == [1, 2]
    @test H2Trees.nodesatlevel(tree) == H2Trees.treeindex(tree).nodes_by_level
    @test H2Trees.treeindex(tree).leaves == H2Trees.leaves(tree)
    @test H2Trees.depthfirstnodes(tree) == collect(Int, H2Trees.DepthFirstIterator(tree))
end

@testset "Float32" begin
    tree = TwoNTree(SVector(0.0f0, 0.0f0, 0.0f0), 1.0f0)

    @test H2Trees.root(tree) == 1
    @test H2Trees.center(tree, 1) == SVector(0.0f0, 0.0f0, 0.0f0)
    @test H2Trees.halfsize(tree, 1) == 1.0f0
    @test H2Trees.level(tree, 1) == 1

    @test H2Trees.treetrait(tree) == H2Trees.isTwoNTree()
end

@testset "TreeIndex invariants" begin
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere4.in")
    )
    points = vertices(m)
    tree = buildtree(points; builder=TwoNTreeBuilder(; minhalfsize=0.1, minvalues=10))

    @test !ismutabletype(H2Trees.TreeIndex)
    @test !ismutabletype(typeof(tree))
    @test fieldtype(typeof(tree), :index) <: Base.RefValue{H2Trees.TreeIndex}

    index = H2Trees.treeindex(tree)
    @test H2Trees.levels(tree) == index.minlevel:index.maxlevel
    for level in H2Trees.levels(tree)
        for node in H2Trees.nodesatlevel(tree, level)
            @test H2Trees.level(tree, node) == level
        end
    end
    @test sort(index.leaves) == sort(H2Trees.leaves(tree))
    @test sort(index.leaves) == sort([
        node for node in H2Trees.DepthFirstIterator(tree) if H2Trees.isleaf(tree, node)
    ])

    # Rebuilding replaces the boxed index with a whole new coherent value, not a mutated one.
    before = H2Trees.treeindex(tree)
    H2Trees.rebuildtreeindex!(tree)
    after = H2Trees.treeindex(tree)
    @test before !== after
    # `TreeIndex` has no custom `==` (its Vector fields make the default `==` fall back to `===`),
    # so compare the coherent value fieldwise instead.
    @test before.nodes_by_level == after.nodes_by_level
    @test before.depthfirstnodes == after.depthfirstnodes
    @test before.leaves == after.leaves
    @test before.minlevel == after.minlevel
    @test before.maxlevel == after.maxlevel

    # Topology-changing operation (balanceleaves! on a BoundingBallTree) leaves the invariants
    # holding after the index is rebuilt.
    balltree = H2Trees.BoundingBallTree(
        [
            H2Trees.Node(
                H2Trees.BoundingBallData(Int[], SVector(0.0, 0.0), 2.0, 1), 0, 0, 2
            ),
            H2Trees.Node(
                H2Trees.BoundingBallData([1], SVector(-1.0, 0.0), 1.0, 2), 3, 1, 0
            ),
            H2Trees.Node(
                H2Trees.BoundingBallData(Int[], SVector(1.0, 0.0), 1.0, 2), 0, 1, 4
            ),
            H2Trees.Node(H2Trees.BoundingBallData([2], SVector(1.0, 0.0), 0.5, 3), 0, 3, 0),
        ],
        1,
        SVector(0.0, 0.0),
        2.0,
        [[1], [2, 3], [4]],
    )
    H2Trees.balanceleaves!(balltree)
    ballindex = H2Trees.treeindex(balltree)
    @test H2Trees.levels(balltree) == ballindex.minlevel:ballindex.maxlevel
    @test sort(ballindex.leaves) == sort(H2Trees.leaves(balltree))
end

@testset "Simple Blocktree" begin
    λ = 1.0
    minhalfsize = λ / 10

    mx = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere5.in")
    )

    my = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "spherewithcenter7.in")
    )

    X = raviartthomas(mx)
    Y = raviartthomas(my)

    tree = buildtree(
        X,
        Y;
        builder=BlockTreeBuilder(;
            test=TwoNTreeBuilder(;
                minhalfsize=minhalfsize, minvalues=10, protrusion=NoProtrusionCheck()
            ),
            trial=TwoNTreeBuilder(;
                minhalfsize=minhalfsize, minvalues=3, protrusion=NoProtrusionCheck()
            ),
        ),
    )

    for tree in [H2Trees.testtree(tree), H2Trees.trialtree(tree)]
        valuesatnodes = H2Trees.valuesatnodes(tree)
        for (functionid, value) in enumerate(valuesatnodes)
            @test length(value) == 1
            @test functionid in H2Trees.values(tree, value[1])
        end
        nodesatvalues = H2Trees.nodesatvalues(tree)
        for (key, value) in nodesatvalues
            @test length(key) == 1
            key = key[1]
            @test sort(value) == sort(H2Trees.values(tree, key))
        end
    end

    @test eltype(tree) == SVector{3,Float64}
    @test eltype(H2Trees.testtree(tree)) == SVector{3,Float64}

    @test H2Trees.treewithmorelevels(tree) == H2Trees.trialtree(tree)

    tree2 = buildtree(
        X,
        Y;
        builder=BlockTreeBuilder(;
            test=TwoNTreeBuilder(;
                minhalfsize=minhalfsize, minvalues=0, protrusion=NoProtrusionCheck()
            ),
            trial=TwoNTreeBuilder(;
                minhalfsize=minhalfsize, minvalues=0, protrusion=NoProtrusionCheck()
            ),
        ),
    )
    @test H2Trees.treewithmorelevels(tree2) == H2Trees.testtree(tree2)

    @test H2Trees.minhalfsize(H2Trees.trialtree(tree)) == minhalfsize

    for level in H2Trees.levels(H2Trees.trialtree(tree))
        halfsize = H2Trees.halfsize(
            H2Trees.trialtree(tree),
            H2Trees.LevelIterator(H2Trees.trialtree(tree), level)[begin],
        )

        for node in H2Trees.LevelIterator(H2Trees.trialtree(tree), level)
            @test H2Trees.level(H2Trees.trialtree(tree), node) == level
            @test H2Trees.halfsize(H2Trees.trialtree(tree), node) == halfsize

            for testnode in H2Trees.LevelIterator(H2Trees.testtree(tree), level)
                @test H2Trees.level(H2Trees.trialtree(tree), node) ==
                    H2Trees.level(H2Trees.testtree(tree), testnode)
                @test H2Trees.halfsize(H2Trees.trialtree(tree), node) ==
                    H2Trees.halfsize(H2Trees.testtree(tree), testnode)
            end
        end

        level == 1 && continue

        level - 1 ∉ H2Trees.levels(H2Trees.trialtree(tree)) && continue

        halfsizeabove = H2Trees.halfsize(
            H2Trees.trialtree(tree),
            H2Trees.LevelIterator(H2Trees.trialtree(tree), level - 1)[begin],
        )

        @test halfsizeabove == 2 * halfsize
    end

    @test H2Trees.treetrait(tree) == H2Trees.isBlockTree()

    display(tree)
    display(H2Trees.testtree(tree))
    display(H2Trees.trialtree(tree))
end

@testset "Compute buffers" begin
    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "sphere6.in")
    )
    points = vertices(m)
    tree = buildtree(
        points; builder=TwoNTreeBuilder(; minhalfsize=0.1, root=2, minlevel=2, minvalues=10)
    )

    buffers = H2Trees.computevectorbuffers(tree, ComplexF64)

    for (node, buffer) in buffers
        @test length(buffer) == length(H2Trees.values(tree, node))
        @test eltype(buffer) == ComplexF64
    end

    my = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "cuboid2.in")
    )

    tree = buildtree(
        raviartthomas(m),
        raviartthomas(my);
        builder=BlockTreeBuilder(;
            test=TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0),
            trial=TwoNTreeBuilder(; minhalfsize=0.1, minvalues=0),
        ),
    )

    testbuffer, trialbuffer = H2Trees.computevectorbuffers(tree, ComplexF64)

    for (node, buffer) in testbuffer
        @test length(buffer) == length(H2Trees.values(H2Trees.testtree(tree), node))
        @test eltype(buffer) == ComplexF64
    end

    for (node, buffer) in trialbuffer
        @test length(buffer) == length(H2Trees.values(H2Trees.trialtree(tree), node))
        @test eltype(buffer) == ComplexF64
    end
end
