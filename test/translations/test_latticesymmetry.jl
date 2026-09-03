using Test
using H2Trees
using Random: Xoshiro
using StaticArrays: SVector
using LinearAlgebra: det

# The far-field offsets that actually occur at one level of a uniform TwoNTree: the children of
# the parent's neighbours, minus the box's own neighbours, unioned over every child slot. The
# union is what matters because `DirectionInvariancePerLevel` deduplicates across a whole level,
# not per box, which is why this is 316 in three dimensions and not the 189 of a single box.
function farfieldoffsets(D::Int)
    slots = Iterators.product(ntuple(_ -> (0, 1), D)...)
    near = Set(Iterators.product(ntuple(_ -> (-1):1, D)...))
    offsets = Set{NTuple{D,Int}}()
    for slot in slots
        for parent in Iterators.product(ntuple(_ -> (-1):1, D)...),
            child in Iterators.product(ntuple(_ -> (0, 1), D)...)

            q = ntuple(i -> 2parent[i] + child[i] - slot[i], D)
            q in near || push!(offsets, q)
        end
    end
    return offsets
end

function orbitcount(offsets, group)
    return length(Set(first(canonicalizetranslation(q, group)) for q in offsets))
end

@testset "Group orders" begin
    @test length(symmetrygroup(FullLatticeSymmetry(), Val(2))) == 8
    @test length(symmetrygroup(FullLatticeSymmetry(), Val(3))) == 48
    @test length(symmetrygroup(AxisPreservingSymmetry(3), Val(3))) == 16
    @test length(symmetrygroup(AxisPreservingSymmetry(1), Val(3))) == 16
    @test length(symmetrygroup(AxisPreservingSymmetry(2), Val(2))) == 4
    @test length(symmetrygroup(OppositeSymmetry(), Val(2))) == 2
    @test length(symmetrygroup(OppositeSymmetry(), Val(3))) == 2
    @test length(symmetrygroup(NoSymmetry(), Val(3))) == 1

    # The identity must be element 1: canonicalization seeds its search with it.
    for D in (2, 3),
        policy in
        (NoSymmetry(), OppositeSymmetry(), AxisPreservingSymmetry(2), FullLatticeSymmetry())

        group = symmetrygroup(policy, Val(D))
        @test group[1] == H2Trees.identitysymmetry(Val(D))
    end

    @test_throws "must lie in 1:2" symmetrygroup(AxisPreservingSymmetry(3), Val(2))
end

@testset "Element algebra" begin
    for D in (2, 3)
        group = symmetrygroup(FullLatticeSymmetry(), Val(D))
        identity = H2Trees.identitysymmetry(Val(D))
        rng = Xoshiro(11)
        vectors = [ntuple(_ -> rand(rng, -4:4), D) for _ in 1:20]

        for s in group
            # Orthogonal: a signed permutation preserves the euclidean norm exactly.
            for q in vectors
                @test sum(abs2, H2Trees.applysymmetry(s, q)) == sum(abs2, q)
            end

            inverse = H2Trees.inversesymmetry(s)
            @test H2Trees.composesymmetry(inverse, s) == identity
            @test H2Trees.composesymmetry(s, inverse) == identity
            for q in vectors
                @test H2Trees.applysymmetry(inverse, H2Trees.applysymmetry(s, q)) == q
            end
        end

        # Composition agrees with successive application, in that order.
        for a in group, b in group
            q = vectors[1]
            @test H2Trees.applysymmetry(H2Trees.composesymmetry(a, b), q) ==
                H2Trees.applysymmetry(a, H2Trees.applysymmetry(b, q))
        end

        # Closure, and the inverse-index table.
        elements = Set(group.elements)
        for a in group, b in group
            @test H2Trees.composesymmetry(a, b) in elements
        end
        for k in eachindex(group)
            @test group[group.inverseindex[k]] == H2Trees.inversesymmetry(group[k])
        end
    end
end

@testset "Constructor validation" begin
    @test_throws "not a permutation" LatticeSymmetry((1, 1, 2), (1, 1, 1))
    @test_throws "outside 1:3" LatticeSymmetry((1, 2, 4), (1, 1, 1))
    @test_throws "must all be +1 or -1" LatticeSymmetry((1, 2, 3), (1, 2, 1))
end

@testset "Canonicalization" begin
    groups = Dict(
        (2, :full) => symmetrygroup(FullLatticeSymmetry(), Val(2)),
        (3, :full) => symmetrygroup(FullLatticeSymmetry(), Val(3)),
        (3, :axis) => symmetrygroup(AxisPreservingSymmetry(3), Val(3)),
        (3, :opposite) => symmetrygroup(OppositeSymmetry(), Val(3)),
        (3, :none) => symmetrygroup(NoSymmetry(), Val(3)),
    )

    @testset "documented examples" begin
        @test first(canonicalizetranslation((-2, 5), groups[(2, :full)])) == (5, 2)
        @test first(canonicalizetranslation((-2, 3, -1), groups[(3, :full)])) == (3, 2, 1)
    end

    @testset "reconstruction invariant" begin
        # q == applysymmetry(group[id], canonical). This is the direction of use and the one a
        # test built only on `OppositeSymmetry` cannot check, the antipodal map being its own
        # inverse, so every group here is exercised, including chiral elements.
        rng = Xoshiro(23)
        for ((D, _), group) in groups, _ in 1:200
            q = ntuple(_ -> rand(rng, -5:5), D)
            canonical, id = canonicalizetranslation(q, group)
            @test H2Trees.applysymmetry(group[id], canonical) == q
        end
    end

    @testset "orbit invariant" begin
        # Every member of an orbit canonicalizes to the same representative. This is the
        # property the deduplication key rests on.
        rng = Xoshiro(29)
        for ((D, _), group) in groups, _ in 1:50
            q = ntuple(_ -> rand(rng, -5:5), D)
            canonical = first(canonicalizetranslation(q, group))
            for r in group
                @test first(canonicalizetranslation(H2Trees.applysymmetry(r, q), group)) ==
                    canonical
            end
        end
    end

    @testset "degenerate orbits are deterministic" begin
        # Offsets on a symmetry axis or plane are fixed by several elements, so the tie-break
        # (lowest-numbered symmetry attaining the representative) is what makes the ID stable.
        degenerate2 = [(0, 0), (3, 0), (0, 3), (3, 3), (-3, 3)]
        degenerate3 = [
            (0, 0, 0), (3, 0, 0), (0, 3, 0), (3, 3, 0), (3, 3, 3), (3, 2, 0), (3, 2, 1)
        ]
        for (D, offsets) in ((2, degenerate2), (3, degenerate3))
            for (_, group) in filter(p -> first(first(p)) == D, groups), q in offsets
                first_call = canonicalizetranslation(q, group)
                @test canonicalizetranslation(q, group) == first_call
                @test H2Trees.applysymmetry(group[first_call[2]], first_call[1]) == q
                @test length(symmetryorbit(q, group)) <= length(group)
            end
        end
    end

    @testset "canonical region of the full group" begin
        # Lexicographic maximization over all signed permutations lands in q1 >= q2 >= ... >= 0.
        rng = Xoshiro(31)
        for D in (2, 3), _ in 1:200
            q = ntuple(_ -> rand(rng, -5:5), D)
            canonical = first(canonicalizetranslation(q, groups[(D, :full)]))
            @test issorted(collect(canonical); rev=true)
            @test canonical[D] >= 0
            @test sort(collect(abs.(q)); rev=true) == collect(canonical)
        end
    end

    @testset "NoSymmetry reproduces plain deduplication" begin
        rng = Xoshiro(37)
        for _ in 1:50
            q = ntuple(_ -> rand(rng, -5:5), 3)
            @test canonicalizetranslation(q, groups[(3, :none)]) == (q, 1)
        end
    end

    @testset "preallocated form agrees" begin
        rng = Xoshiro(41)
        for ((D, _), group) in groups, _ in 1:100
            q = ntuple(_ -> rand(rng, -5:5), D)
            canonical = Vector{Int}(undef, D)
            scratch = Vector{Int}(undef, D)
            id = canonicalizetranslation!(canonical, scratch, collect(q), group)
            @test (NTuple{D,Int}(canonical), id) == canonicalizetranslation(q, group)
        end
    end
end

@testset "A larger group's representative is not reachable in a smaller one" begin
    # The reason canonicalization must run over the subgroup the consumer can actually act with.
    # An element preserving axis 3 can only draw the third component from the third component, so
    # it fixes |q[3]|. For (1, 2, 3) the full group's representative (3, 2, 1) has a third
    # component of 1 and is therefore unreachable: a consumer restricted to that subgroup could
    # not transform it back into the offset it needs.
    #
    # Note this is NOT visible for every offset. (-2, 3, -1) also canonicalizes to (3, 2, 1)
    # under the full group and IS reachable, because its polar component already carries the
    # smallest magnitude. The mismatch is silent on exactly the offsets that happen to be
    # aligned, which is what makes it worth pinning.
    full = symmetrygroup(FullLatticeSymmetry(), Val(3))
    axis = symmetrygroup(AxisPreservingSymmetry(3), Val(3))

    q = (1, 2, 3)
    fullcanonical = first(canonicalizetranslation(q, full))
    @test fullcanonical == (3, 2, 1)
    @test !any(H2Trees.applysymmetry(s, fullcanonical) == q for s in axis)
    @test any(H2Trees.applysymmetry(s, fullcanonical) == (-2, 3, -1) for s in axis)

    # Canonicalizing over the subgroup instead gives a representative that is reachable.
    axiscanonical, id = canonicalizetranslation(q, axis)
    @test H2Trees.applysymmetry(axis[id], axiscanonical) == q
end

@testset "Translation reduction on the far-field offsets" begin
    # The payoff, as a regression rather than a claim. These counts are the reason the
    # three-dimensional feature is scoped around the order-16 subgroup: the full group would be
    # worth 19.75x but a tensor-product cos(theta)-by-phi grid cannot act with it, and the
    # antipodal map alone, the obvious starting point, is worth only 2x.
    offsets2 = farfieldoffsets(2)
    offsets3 = farfieldoffsets(3)
    @test length(offsets2) == 40
    @test length(offsets3) == 316

    @test orbitcount(offsets2, symmetrygroup(FullLatticeSymmetry(), Val(2))) == 7
    @test orbitcount(offsets2, symmetrygroup(OppositeSymmetry(), Val(2))) == 20

    @test orbitcount(offsets3, symmetrygroup(FullLatticeSymmetry(), Val(3))) == 16
    @test orbitcount(offsets3, symmetrygroup(AxisPreservingSymmetry(3), Val(3))) == 34
    @test orbitcount(offsets3, symmetrygroup(OppositeSymmetry(), Val(3))) == 158
    @test orbitcount(offsets3, symmetrygroup(NoSymmetry(), Val(3))) == 316

    # Any policy is a subgroup of the full one, so no policy can beat it.
    for D in (2, 3)
        offsets = D == 2 ? offsets2 : offsets3
        full = orbitcount(offsets, symmetrygroup(FullLatticeSymmetry(), Val(D)))
        for policy in (NoSymmetry(), OppositeSymmetry(), AxisPreservingSymmetry(D))
            @test orbitcount(offsets, symmetrygroup(policy, Val(D))) >= full
        end
    end
end

@testset "Symmetry-reduced translation trait" begin
    # The three properties the trait must have on a real tree: the interaction graph is
    # UNCHANGED, the unique translation count drops, and every interaction can recover its own
    # direction from the canonical one it points at.
    rng = Xoshiro(7)
    points = [SVector{3,Float64}(randn(rng), randn(rng), randn(rng)) for _ in 1:4000]
    tree = H2Trees.TwoNTree(
        points; builder=H2Trees.TwoNTreeBuilder(; minhalfsize=0.25, minvalues=0)
    )
    @test H2Trees.numberoflevels(tree) > 3

    plans = H2Trees.buildplans(tree)
    plan = H2Trees.translatingplan(plans.trialaggregationplan, plans.testdisaggregationplan)

    plaininfos, plaindirections, _ = H2Trees.translations(
        tree, plan, H2Trees.DirectionInvariancePerLevel()
    )
    plainpairs = sum(length, plaininfos)

    # The default policy is the conservative one.
    @test SymmetryDirectionInvariancePerLevel().policy === OppositeSymmetry()

    for policy in (OppositeSymmetry(), AxisPreservingSymmetry(3), FullLatticeSymmetry())
        infos, directions, _ = H2Trees.translations(
            tree, plan, SymmetryDirectionInvariancePerLevel(policy)
        )
        group = symmetrygroup(policy, Val(3))

        # Same interaction graph: symmetry changes what is STORED, never what is computed.
        @test sum(length, infos) == plainpairs
        # Never worse than plain deduplication, and strictly better than the trivial group.
        @test length(directions) <= length(plaindirections)
        policy === OppositeSymmetry() && @test length(directions) < length(plaindirections)

        # The reconstruction invariant, on every interaction rather than a sample.
        allrecover = true
        for leveltfinfo in infos, tfinfo in leveltfinfo
            h = H2Trees.halfsize(tree, tfinfo.receivingnode)
            d =
                H2Trees.center(tree, tfinfo.receivingnode) -
                H2Trees.center(tree, tfinfo.translatingnode)
            actual = ntuple(i -> round(Int, d[i] / h), 3)
            canonical = ntuple(i -> round(Int, directions[tfinfo.translationID][i] / h), 3)
            if H2Trees.applysymmetry(group[tfinfo.symmetryID], canonical) != actual
                allrecover = false
                break
            end
        end
        @test allrecover
    end

    # Antipodal reduction is exactly 2x on a tree with no self-opposite offsets (q = -q only for
    # q = 0, which is never a translation), so this is a hard number rather than a bound.
    _, antipodal, _ = H2Trees.translations(
        tree, plan, SymmetryDirectionInvariancePerLevel(OppositeSymmetry())
    )
    @test 2 * length(antipodal) == length(plaindirections)
end

# TWO TREES, WHERE THE OFFSET DECIDES WHETHER A SYMMETRY IS APPLICABLE AT ALL.
#
# A block tree's roots can sit anywhere relative to each other, and the plain
# `DirectionInvariancePerLevel` handles that by folding the root-center offset out before taking the
# integer key: every interaction sharing a key then has the same physical displacement
# `halfsize * key + offset`, so one stored translation serves them all.
#
# A SYMMETRY MERGES DIFFERENT KEYS, and that is where the offset stops being harmless. The stored
# direction is `halfsize * canonical + offset` and the interaction needs
# `halfsize * (S * canonical) + offset`, but applying `S` yields `halfsize * (S * canonical) +
# S * offset`. Those agree only when `S * offset == offset`, which for a generic offset means the
# identity alone, so the reduction would be wrong, and quietly wrong.
#
# The implemented condition is therefore that the displacements lie on a lattice, which admits
# exactly the two correct cases (offset zero, offset a lattice vector) and refuses the third. All
# three are pinned below, because the failure has no symptom other than a slightly wrong field.
@testset "block trees: the root offset decides applicability" begin
    rng = Xoshiro(0xABCDEF)
    cloud = [SVector{3,Float64}(rand(rng, 3)) for _ in 1:3000]
    builder = TwoNTreeBuilder(; minhalfsize=0.1, minvalues=1)
    policy = AxisPreservingSymmetry(3)
    group = H2Trees.symmetrygroup(policy, Val(3))
    identity3 = H2Trees.identitysymmetry(Val(3))

    roothalfsize = H2Trees.halfsize(
        TwoNTree(cloud; builder=builder), H2Trees.root(TwoNTree(cloud; builder=builder))
    )

    blocktree(shift) = H2Trees.buildtree(
        [p + SVector{3,Float64}(shift, 0.0, 0.0) for p in cloud],
        cloud;
        builder=BlockTreeBuilder(; test=builder, trial=builder),
    )

    function translatingplanof(tree)
        plans = buildplans(tree; builder=PlanBuilder())
        return H2Trees.translatingplan(
            plans.trialaggregationplan, plans.testdisaggregationplan
        )
    end

    # Reconstruction on EVERY interaction, computed from the two trees' own centers rather than
    # from anything `translations` returned besides the tag it is being checked against.
    function allreconstruct(tree, infos, directions)
        testtree = H2Trees.testtree(tree)
        trialtree = H2Trees.trialtree(tree)
        for leveltfinfo in infos, tfinfo in leveltfinfo
            h = H2Trees.halfsize(testtree, tfinfo.receivingnode)
            d =
                H2Trees.center(testtree, tfinfo.receivingnode) -
                H2Trees.center(trialtree, tfinfo.translatingnode)
            actual = ntuple(i -> round(Int, d[i] / h), 3)
            canonical = ntuple(i -> round(Int, directions[tfinfo.translationID][i] / h), 3)
            H2Trees.applysymmetry(group[tfinfo.symmetryID], canonical) == actual ||
                return false
        end
        return true
    end

    nonidentitycount(infos) = sum(
        sum(1 for i in lvl if group[i.symmetryID] != identity3; init=0) for lvl in infos
    )

    @testset "coincident roots (offset zero)" begin
        tree = H2Trees.buildtree(
            cloud, cloud; builder=BlockTreeBuilder(; test=builder, trial=builder)
        )
        plan = translatingplanof(tree)
        _, plaindirections, _ = H2Trees.translations(
            tree, plan, H2Trees.DirectionInvariancePerLevel()
        )
        infos, directions, _ = H2Trees.translations(
            tree, plan, SymmetryDirectionInvariancePerLevel(policy)
        )
        @test length(directions) < length(plaindirections)
        @test nonidentitycount(infos) > 0
        @test allreconstruct(tree, infos, directions)
    end

    @testset "roots offset by a lattice vector" begin
        # One root halfsize. Finer halfsizes divide it, so the displacement stays a lattice vector
        # at every level, which is the property the check needs and is not implied by alignment at the
        # finest level alone, since the halfsize doubles going up.
        #
        # A larger shift separates the clouds entirely and leaves a single far interaction, which
        # is how the Petrov testset in MLFMA came to assert reconstruction on a case where no
        # symmetry ever ran. The counts below are what keeps this one honest.
        tree = blocktree(roothalfsize)
        plan = translatingplanof(tree)
        _, plaindirections, _ = H2Trees.translations(
            tree, plan, H2Trees.DirectionInvariancePerLevel()
        )
        infos, directions, _ = H2Trees.translations(
            tree, plan, SymmetryDirectionInvariancePerLevel(policy)
        )
        @test sum(length, infos) > 1000                 # not the one-interaction trap
        @test length(directions) < length(plaindirections) / 4
        @test nonidentitycount(infos) > 1000            # symmetries genuinely applied
        @test allreconstruct(tree, infos, directions)
    end

    @testset "roots offset by a non-lattice vector are refused" begin
        tree = blocktree(4 * roothalfsize + 0.37 * roothalfsize)
        plan = translatingplanof(tree)
        # The plain trait still works: it folds the offset out and never merges keys, so it is
        # correct for any offset. That contrast is the point: this is a restriction on symmetry,
        # not on block trees.
        _, plaindirections, _ = H2Trees.translations(
            tree, plan, H2Trees.DirectionInvariancePerLevel()
        )
        @test !isempty(plaindirections)

        # A DISTINCT EXCEPTION TYPE, not just a message. Consumers are meant to catch this and
        # fall back to the plain trait; a broad `catch` there would swallow method errors and
        # out-of-memory as "not lattice aligned" and hand back an unreduced collection with no
        # indication why. Both the type and the text are pinned: the type is the contract, the
        # text is what a user without a fallback actually reads.
        @test_throws H2Trees.NonLatticeTranslationError H2Trees.translations(
            tree, plan, SymmetryDirectionInvariancePerLevel(policy)
        )
        @test_throws "not a multiple of the level halfsize" H2Trees.translations(
            tree, plan, SymmetryDirectionInvariancePerLevel(policy)
        )
    end
end

@testset "determinant: rotation vs reflection" begin
    # A signed coordinate permutation is an orthogonal matrix, so its determinant is +/-1. The value
    # matters to a consumer because a pseudo-tensor coupling (one built from a curl or cross
    # product, such as a double-layer kernel) transforms as `det(S) · S K Sᵀ` where a true tensor
    # transforms as `S K Sᵀ`. The two are indistinguishable by shape, and differ by exactly a sign.
    for D in (2, 3)
        group = H2Trees.symmetrygroup(H2Trees.FullLatticeSymmetry(), Val(D))

        # Against the actual matrix, so this tests the integer shortcut rather than restating it.
        for symmetry in group
            S = zeros(Float64, D, D)
            for i in 1:D
                S[i, symmetry.permutation[i]] = symmetry.signs[i]
            end
            @test H2Trees.determinant(symmetry) == round(Int, det(S))
        end

        # HALF the group are reflections, which is why this cannot be discovered by testing the
        # identity or any pure rotation: a consumer applying the wrong rule is exactly right on one
        # half and O(1) wrong on the other.
        determinants = [H2Trees.determinant(s) for s in group]
        @test count(==(1), determinants) == length(group) ÷ 2
        @test count(==(-1), determinants) == length(group) ÷ 2

        # It is a group homomorphism, which is what makes it safe to apply per element.
        for a in group, b in group
            @test H2Trees.determinant(H2Trees.composesymmetry(a, b)) ==
                H2Trees.determinant(a) * H2Trees.determinant(b)
        end
        @test H2Trees.determinant(H2Trees.identitysymmetry(Val(D))) == 1
        for symmetry in group
            @test H2Trees.determinant(H2Trees.inversesymmetry(symmetry)) ==
                H2Trees.determinant(symmetry)
        end
    end
end

# THE GEOMETRY CONTRACT, AND WHY IT IS A DIFFERENT KIND OF REFUSAL FROM THE OFFSET ONE ABOVE.
#
# `SymmetryDirectionInvariancePerLevel` needs a LATTICE: a symmetry maps a displacement to another
# displacement on the same lattice, which requires centers separated by integer multiples of the
# level halfsize. That is a property of the regular 2^D subdivision, not of trees in general.
#
# Two failures are possible and they are deliberately NOT the same type:
#
#   * wrong geometry, an `ArgumentError`. A configuration error: no offset makes a ball
#     tree lattice-aligned, so there is nothing to fall back from and nothing should catch it.
#   * right geometry but the wrong offset, a `NonLatticeTranslationError`. A fallback condition, which
#     consumers ARE expected to catch and answer with `DirectionInvariancePerLevel` (MLFMA does).
#
# Before the explicit methods this testset pins, the first case was a `MethodError` naming
# `_translations` and a list of functors: the generic `translations(tree, ::AbstractTreeTrait, ...)`
# forwards seven positional arguments while the symmetry trait's `_translations` takes eight. That
# is a dispatch gap, not a contract, and it is indistinguishable from a genuine bug at the call site.
@testset "the symmetry trait refuses non-lattice GEOMETRY, distinctly from a non-lattice OFFSET" begin
    # A deterministic median bisection along the widest axis. Written here rather than reaching for
    # `KMeansTree`, which lives behind an extension: all this fixture has to be is a real, plan-able
    # tree that is not a `TwoNTree`.
    function bisect(points, globalpointids, level, numsplits; kwargs...)
        pts = [points[i] for i in globalpointids]
        axis = argmax([
            maximum(p[d] for p in pts) - minimum(p[d] for p in pts) for d in 1:3
        ])
        ids = globalpointids[sortperm(globalpointids; by=i -> points[i][axis])]
        half = max(1, length(ids) ÷ 2)
        groups = filter(!isempty, [ids[1:half], ids[(half + 1):end]])
        centers = [sum(points[i] for i in g) / length(g) for g in groups]
        radii = [
            maximum(sqrt(sum(abs2, points[i] - c)) for i in g) for
            (g, c) in zip(groups, centers)
        ]
        return groups, centers, radii
    end

    rng = Xoshiro(0xBEEF)
    cloud = [SVector{3,Float64}(rand(rng, 3)) for _ in 1:2000]
    balltree = H2Trees.buildtree(
        cloud;
        builder=H2Trees.BoundingBallTreeBuilder(;
            splitter=bisect, numsplits=2, minvalues=20
        ),
    )
    twontree = H2Trees.buildtree(
        cloud; builder=TwoNTreeBuilder(; minhalfsize=0.1, minvalues=20)
    )
    @test H2Trees.treetrait(balltree) === H2Trees.isBoundingBallTree()

    ballplan = let p = H2Trees.buildplans(balltree)
        H2Trees.translatingplan(p.trialaggregationplan, p.testdisaggregationplan)
    end
    twonplan = let p = H2Trees.buildplans(twontree)
        H2Trees.translatingplan(p.trialaggregationplan, p.testdisaggregationplan)
    end
    sym = SymmetryDirectionInvariancePerLevel(OppositeSymmetry())

    # NON-VACUITY: the ball tree really does support the plain trait, so what follows is the
    # SYMMETRY being refused and not the fixture being unusable.
    _, balldirections, _ = H2Trees.translations(
        balltree, ballplan, H2Trees.DirectionInvariancePerLevel()
    )
    @test !isempty(balldirections)

    err = try
        H2Trees.translations(balltree, ballplan, sym)
        nothing
    catch exception
        exception
    end
    @test err isa ArgumentError
    # The two refusals must stay distinguishable: a consumer catching the fallback condition must
    # NOT swallow the configuration error.
    @test !(err isa H2Trees.NonLatticeTranslationError)
    @test occursin("TwoNTree", err.msg)
    @test occursin("DirectionInvariancePerLevel", err.msg)

    # ...and the supported geometry is untouched by the new method.
    _, symdirections, _ = H2Trees.translations(twontree, twonplan, sym)
    @test !isempty(symdirections)

    # THE BLOCK-TREE TWIN, called through the trait-dispatched form ON PURPOSE. `BlockTreeBuilder`
    # already rejects a non-`TwoNTree` builder (`_validate_block_tree_builders` accepts only a pair
    # of `TwoNTreeBuilder`s) and `BlockTree(a, b)` builds `TwoNTree`s from POSITIONS, so a mixed
    # block tree cannot be constructed through any public path today. The refusal is therefore
    # defence-in-depth behind that validation rather than a live guard, reachable only by
    # dispatching the traits directly, which is exactly what this does, or by a future tree type
    # that becomes block-able.
    for (testtrait, trialtrait, a, b) in (
        (H2Trees.isBoundingBallTree(), H2Trees.isTwoNTree(), balltree, twontree),
        (H2Trees.isTwoNTree(), H2Trees.isBoundingBallTree(), twontree, balltree),
        (H2Trees.isBoundingBallTree(), H2Trees.isBoundingBallTree(), balltree, balltree),
    )
        blockerr = try
            H2Trees.translations(a, b, testtrait, trialtrait, twonplan, sym)
            nothing
        catch exception
            exception
        end
        @test blockerr isa ArgumentError
        @test !(blockerr isa H2Trees.NonLatticeTranslationError)
    end

    # The all-`TwoNTree` pair must still reach the real implementation: the refusal is a fallback,
    # so a method that accidentally shadowed the specific one would disable the feature entirely.
    _, blockdirections, _ = H2Trees.translations(
        twontree, twontree, H2Trees.isTwoNTree(), H2Trees.isTwoNTree(), twonplan, sym
    )
    @test !isempty(blockdirections)
end

# THE ERROR CARRIES THE TREE'S OWN COORDINATE TYPE.
#
# `NonLatticeTranslationError`'s fields used to be `Vector{Float64}`/`Float64`. A review read that as
# breaking the catchable-type contract for non-`Float64` trees, giving a `MethodError` in place of the
# intended exception. IT DID NOT: Julia's default constructor converts each field, so `Float32`,
# `BigFloat` and `Rational` all produced a genuine, catchable `NonLatticeTranslationError`. The first
# two assertions below pin that, because it is the property consumers actually depend on and nothing
# else states it.
#
# WHAT THE WIDENING DID BREAK IS THE MESSAGE, which is what the rest of this testset is about. A
# `Float32` tree whose halfsize is `0.8f0` reported it as `0.800000011920929`, the binary expansion
# of a number the user never wrote, and a `BigFloat` tree lost precision outright. A type not
# convertible to `Float64` at all would have raised a conversion error instead of this one.
#
# THE 20 ASSERTIONS BELOW FALL INTO THREE KINDS, and mixing them up is how a testset comes to look
# more rigorous than it is. Measured by reverting the struct and rerunning, not asserted from
# reading:
#
#   1. CONTRACT (5, pass before AND after). `err isa NonLatticeTranslationError`,
#      `!(err isa ArgumentError)`, and the `e isa ...` check inside the three-type loop. These pin
#      the catchable-error contract, the thing consumers depend on and the thing a review
#      doubted, but they are not regression tests for the parameterization, because the pre-change
#      struct satisfied them too.
#   2. REGRESSION (9, fail before). `err.halfsize isa Float32`, `eltype(err.translation) === Float32`,
#      the message not containing `0.800000011920929`, and the six struct-level preservation checks
#      across `Float32`/`BigFloat`/`Rational`. These are what actually separate the two versions.
#   3. VALUE (6, pass before AND after). `roothalfsize isa Float32`, `err.halfsize == roothalfsize`,
#      the message containing `0.8`, and the three `showerror` smoke checks. Useful, since they
#      guard against carrying the wrong number or a broken message, but blind to the storage type,
#      because `==` promotes and `Float64(0.8f0)` is exactly the widened value.
#
# 5 + 9 + 6 = 20, and the 9 is the count the revert actually produced. The loop contributes FOUR
# assertions per type, not three, which is what made two earlier drafts of this comment miscount.
#
# The moral for anyone extending this: the obvious "build a Float32 tree and check the error is
# catchable" test is category 1. It passes either way and would have asserted nothing about the
# change.
@testset "the non-lattice error keeps the tree's coordinate type" begin
    rng = Xoshiro(0xABCDEF)
    cloud = [SVector{3,Float32}(rand(rng, Float32, 3)) for _ in 1:3000]
    builder = TwoNTreeBuilder(; minhalfsize=0.1f0, minvalues=1)
    roothalfsize = let t = TwoNTree(cloud; builder=builder)
        H2Trees.halfsize(t, H2Trees.root(t))
    end
    @test roothalfsize isa Float32

    # Deliberately not a lattice vector, matching the offset case pinned above; the fractional
    # part is what the check rejects.
    shift = 4roothalfsize + 0.37f0 * roothalfsize
    tree = H2Trees.buildtree(
        [p + SVector{3,Float32}(shift, 0, 0) for p in cloud],
        cloud;
        builder=BlockTreeBuilder(; test=builder, trial=builder),
    )
    plans = buildplans(tree; builder=PlanBuilder())
    plan = H2Trees.translatingplan(plans.trialaggregationplan, plans.testdisaggregationplan)

    err = try
        H2Trees.translations(
            tree, plan, SymmetryDirectionInvariancePerLevel(OppositeSymmetry())
        )
        nothing
    catch exception
        exception
    end

    # The contract the review doubted, which held all along and now has a test.
    @test err isa H2Trees.NonLatticeTranslationError
    @test !(err isa ArgumentError)

    # The part that genuinely was wrong: the fields are the TREE's type, not `Float64`.
    @test err.halfsize isa Float32
    @test eltype(err.translation) === Float32
    # The VALUE is the tree's own halfsize, not merely something of the right type. NOTE this one
    # does NOT discriminate between the two struct versions and is not meant to: `==` promotes, and
    # `Float64(0.8f0)` is exactly the `0.800000011920929` the widened field used to hold, so it
    # passed before the parameterization too. It guards a different regression: a future change
    # that preserved the type while carrying the wrong number.
    @test err.halfsize == roothalfsize

    # ...and therefore the message shows the user their own number. `0.800000011920929` is what the
    # widened field printed, and seeing it is how this was found.
    msg = sprint(showerror, err)
    @test occursin("0.8", msg)
    @test !occursin("0.800000011920929", msg)

    # Type preservation is a property of the struct, not of one tree that happens to be Float32.
    for (T, tr, hs) in (
        (Float32, Float32[1.5, 2.5, 3.5], 0.5f0),
        (BigFloat, BigFloat[1.5, 2.5], BigFloat(0.5)),
        (Rational{Int}, Rational{Int}[3 // 2, 5 // 2], 1//2),
    )
        e = H2Trees.NonLatticeTranslationError(7, 9, tr, hs, 1, zero(hs))
        @test e isa H2Trees.NonLatticeTranslationError
        @test e.halfsize isa T
        @test eltype(e.translation) === T
        # `showerror` divides `translation[coordinate]` by `halfsize`; for `Rational` that stays
        # exact, which is the case a `Float64` field could not have represented at all.
        @test occursin("halfsize", sprint(showerror, e))
    end
end
