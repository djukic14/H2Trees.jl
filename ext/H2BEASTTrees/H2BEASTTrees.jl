module H2BEASTTrees
using StaticArrays
using BEAST
using H2Trees
using Graphs

import H2Trees: TwoNTree, BlockTree, isgalerkinsymmetric
import H2Trees: BEASTProtrusionFunctor, halfsize, center, adjacencygraph
import H2Trees: defaultprotrusioncheck, defaultprotrusionfunctor
import H2Trees: KMeansTree, MetisForest, MetisTree
import H2Trees: buildtree, buildforest, TwoNTreeBuilder, BlockTreeBuilder
import H2Trees: KMeansTreeBuilder, MetisTreeBuilder, MetisForestBuilder

"""
    defaultprotrusionfunctor(space::BEAST.Space)

Return the BEAST-aware protrusion functor for `space`.

This is used when a [`ProtrusionCheck`](@ref) carries the generic
`ComputeProtrusionFunctor()` sentinel: BEAST spaces replace that generic point
functor with `BEASTProtrusionFunctor(space)` so construction measures protrusion
against the actual basis-function support.
"""
defaultprotrusionfunctor(space::BEAST.Space) = BEASTProtrusionFunctor(space)

"""
    defaultprotrusioncheck(space::BEAST.Space)

Return the default protrusion policy for a BEAST space.

`TwoNTreeBuilder()` and each side of `BlockTreeBuilder()` default to
`AutoProtrusionCheck()`. For a `BEAST.Space`, that resolves to
`ProtrusionCheck(; max=0.25, compute=BEASTProtrusionFunctor(space))`, so
`TwoNTree(space)`/`buildtree(space)` work out of the box with the BEAST
protrusion check enabled. Pass `NoProtrusionCheck()` explicitly to opt out, or
pass a `ProtrusionCheck` explicitly to choose a different maximum.
"""
function defaultprotrusioncheck(space::BEAST.Space)
    return H2Trees.ProtrusionCheck(; max=0.25, compute=BEASTProtrusionFunctor(space))
end

"""
    TwoNTree(space::BEAST.Space; builder=TwoNTreeBuilder())

Build a `TwoNTree` from the function positions of a BEAST space.

The default builder uses the BEAST protrusion policy from
[`defaultprotrusioncheck`](@ref). Passing `NoProtrusionCheck()` opts out, and an
explicit `ProtrusionCheck` keeps its `max` while replacing only the generic
point protrusion functor with `BEASTProtrusionFunctor(space)`.
"""
function TwoNTree(space::BEAST.Space; builder=TwoNTreeBuilder())
    return buildtree(space; builder=builder)
end

"""
    KMeansTree(space::BEAST.Space; builder=KMeansTreeBuilder())

Build a `KMeansTree` from `BEAST.positions(space)`.

K-means trees use the point positions directly and do not apply a protrusion check.
"""
function KMeansTree(space::BEAST.Space; builder=KMeansTreeBuilder())
    return buildtree(space; builder=builder)
end

"""
    BlockTree(testspace::BEAST.Space, trialspace::BEAST.Space; builder=BlockTreeBuilder())

Build a Petrov/block tree from a pair of BEAST spaces.

Each side resolves `AutoProtrusionCheck()` against its own BEAST space before the
corresponding `TwoNTree` is built.
"""
function BlockTree(
    testspace::BEAST.Space, trialspace::BEAST.Space; builder=BlockTreeBuilder()
)
    return buildtree(testspace, trialspace; builder=builder)
end

"""
    buildtree(space::BEAST.Space; builder=TwoNTreeBuilder(), graphweights=nothing)

Build a tree from a BEAST space using an explicit builder.

`TwoNTreeBuilder` resolves BEAST protrusion defaults as described by
[`defaultprotrusioncheck`](@ref). `KMeansTreeBuilder` uses only
`BEAST.positions(space)`.
`MetisTreeBuilder` is accepted through this extension's
METIS support; `graphweights` is only valid with that builder.
"""
function buildtree(space::BEAST.Space; builder=TwoNTreeBuilder(), graphweights=nothing)
    return _buildbeasttree(space, builder, graphweights)
end

_checkunusedgraphweights(::Nothing) = nothing
function _checkunusedgraphweights(graphweights)
    return throw(ArgumentError("graphweights is only supported with MetisTreeBuilder"))
end

function _buildbeasttree(space::BEAST.Space, builder::TwoNTreeBuilder, graphweights)
    _checkunusedgraphweights(graphweights)
    builder = H2Trees._resolve_builder_protrusion(builder, space)
    return buildtree(BEAST.positions(space); builder=builder)
end

function _buildbeasttree(space::BEAST.Space, builder::KMeansTreeBuilder, graphweights)
    _checkunusedgraphweights(graphweights)
    return buildtree(BEAST.positions(space); builder=builder)
end

function _buildbeasttree(space::BEAST.Space, builder, graphweights)
    return throw(
        ArgumentError(
            "buildtree does not know how to build from a BEAST.Space with a " *
            "$(typeof(builder)); pass a TwoNTreeBuilder, KMeansTreeBuilder, or MetisTreeBuilder",
        ),
    )
end

"""
    buildtree(testspace::BEAST.Space, trialspace::BEAST.Space; builder=BlockTreeBuilder())

Build a block tree from a BEAST test/trial space pair using `builder`.

The test and trial builders resolve `AutoProtrusionCheck()` independently from
their corresponding BEAST spaces before delegating to the point-based
builder workflow.
"""
function buildtree(
    testspace::BEAST.Space, trialspace::BEAST.Space; builder=BlockTreeBuilder()
)
    builder = BlockTreeBuilder(;
        test=H2Trees._resolve_builder_protrusion(builder.test, testspace),
        trial=H2Trees._resolve_builder_protrusion(builder.trial, trialspace),
    )
    return buildtree(
        BEAST.positions(testspace), BEAST.positions(trialspace); builder=builder
    )
end

include("protrusion.jl")

# Helmholtz3D ##############################################################################

function isgalerkinsymmetric(::Type{<:BEAST.HH3DSingleLayerFDBIO})
    return true
end

function isgalerkinsymmetric(::Type{<:BEAST.HH3DHyperSingularFDBIO})
    return true
end

# Maxwell3D ################################################################################

function isgalerkinsymmetric(::Type{<:BEAST.MWSingleLayer3D})
    return true
end

function isgalerkinsymmetric(::Type{<:BEAST.MWDoubleLayer3D})
    return true
end

# Metis ####################################################################################
include("metis.jl")
end # module H2BEASTTrees
