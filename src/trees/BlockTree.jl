"""
    BlockTree{T}

Container holding a pair of trees used as test and trial clusters.

`BlockTree` is the tree shape used for Petrov-Galerkin plans: the test and
trial sides are built separately but adjusted to a compatible level scale.
"""
struct BlockTree{T}
    testcluster::T
    trialcluster::T
end

"""
    BlockTree(testpositions, trialpositions; builder::BlockTreeBuilder=BlockTreeBuilder())

Construct a `BlockTree` from separate test and trial positions.

Prefer the canonical [`buildtree`](@ref) entry point; this method is the
builder-backed constructor it forwards to. Side-specific protrusion defaults are
resolved before the two inner [`TwoNTree`](@ref)s are built.
"""
function BlockTree(
    testpositions, trialpositions; builder::BlockTreeBuilder=BlockTreeBuilder()
)
    builder = BlockTreeBuilder(;
        test=_resolve_builder_protrusion(builder.test, testpositions),
        trial=_resolve_builder_protrusion(builder.trial, trialpositions),
    )
    return _blocktwontree(testpositions, trialpositions, builder)
end

"""
    _blocktwontree(testpositions, trialpositions, builder::BlockTreeBuilder)

Build a [`BlockTree`](@ref) of two [`TwoNTree`](@ref)s.

Each side is built with compatible root sizes and minimum levels so both trees
can be used together in block-tree traversals.
"""
function _blocktwontree(testpositions, trialpositions, builder::BlockTreeBuilder)
    testcenter, testhalfsize = boundingbox(testpositions)
    trialcenter, trialhalfsize = boundingbox(trialpositions)

    # Promote an untyped `minhalfsize` default (e.g. the `0` from `TwoNTreeBuilder()`) to the
    # coordinate type so it matches the inner `TwoNTree` constructor.
    minhalfsize = oftype(testhalfsize, builder.test.minhalfsize)
    _, testroothalfsize, testminlevel, trialroothalfsize, trialminlevel, = adjusttwontreeblocktreeparameters(
        testhalfsize, trialhalfsize, minhalfsize
    )

    _check_resolved_block_minlevel(builder.test.minlevel, testminlevel, :test)
    _check_resolved_block_minlevel(builder.trial.minlevel, trialminlevel, :trial)

    testtree = TwoNTree(
        SVector(testcenter...),
        testpositions,
        testroothalfsize,
        minhalfsize;
        minlevel=testminlevel,
        root=builder.test.root,
        minvalues=builder.test.minvalues,
        protrusion=builder.test.protrusion,
    )

    trialtree = TwoNTree(
        SVector(trialcenter...),
        trialpositions,
        trialroothalfsize,
        minhalfsize;
        minlevel=trialminlevel,
        root=builder.trial.root,
        minvalues=builder.trial.minvalues,
        protrusion=builder.trial.protrusion,
    )

    # Each `TwoNTree(...)` call above already rebuilds its own tree index once as part of
    # construction: rebuilding again here would just redo that same O(n log n) traversal
    # against an unchanged tree, on both sides.
    return BlockTree(testtree, trialtree)
end

"""
    _check_resolved_block_minlevel(requested, derived, side)

Validate an explicitly requested side minlevel after block-tree geometry has
resolved the required value.
"""
function _check_resolved_block_minlevel(::AutoMinLevel, derived::Int, side::Symbol)
    return nothing
end

function _check_resolved_block_minlevel(requested::Int, derived::Int, side::Symbol)
    requested == derived && return nothing
    return throw(
        ArgumentError(
            "BlockTreeBuilder $side minlevel resolves to $derived from the block-tree level scale, got $requested",
        ),
    )
end

"""
    adjusttwontreeblocktreeparameters(testhalfsize, trialhalfsize, minhalfsize)

Return compatible root halfsizes and minimum levels for the two sides of a
`BlockTree`.

The smaller side may start at a deeper minimum level so both trees share the
same geometric level scale where interactions are compared.
"""
function adjusttwontreeblocktreeparameters(testhalfsize, trialhalfsize, minhalfsize)
    if abs(testhalfsize - trialhalfsize) / max(testhalfsize, trialhalfsize) <
        eps(minhalfsize) * 1e6
        # case where both are the same size
        commonhalfsize = max(testhalfsize, trialhalfsize)

        return minhalfsize,
        roothalfsize(commonhalfsize, minhalfsize),
        1,
        roothalfsize(commonhalfsize, minhalfsize),
        1

    elseif trialhalfsize > testhalfsize
        # case where test tree is larger than trial tree
        testroothalfsize = roothalfsize(testhalfsize, minhalfsize)
        trialroothalfsize = roothalfsize(trialhalfsize, testroothalfsize)

        testminlevel = numberoflevels(trialroothalfsize, testroothalfsize) + 1
        trialminlevel = 1

        return minhalfsize, testroothalfsize, testminlevel, trialroothalfsize, trialminlevel
    else
        # case where trial tree is larger than test tree
        trialroothalfsize = roothalfsize(trialhalfsize, minhalfsize)
        testroothalfsize = roothalfsize(testhalfsize, trialroothalfsize)

        trialminlevel = numberoflevels(testroothalfsize, trialroothalfsize) + 1
        testminlevel = 1

        return minhalfsize, testroothalfsize, testminlevel, trialroothalfsize, trialminlevel
    end
end

"""
    treetrait(::Type{<:BlockTree})

Return [`isBlockTree`](@ref).
"""
function treetrait(::Type{BlockTree{T}}) where {T}
    return isBlockTree()
end

"""
    testtree(tree::BlockTree)

Return the test-side tree.
"""
function testtree(tree::BlockTree)
    return tree.testcluster
end

"""
    trialtree(tree::BlockTree)

Return the trial-side tree.
"""
function trialtree(tree::BlockTree)
    return tree.trialcluster
end

"""
    eltype(tree::BlockTree)

Return the element type of the test-side tree.
"""
function Base.eltype(tree::BlockTree)
    return Base.eltype(testtree(tree))
end
