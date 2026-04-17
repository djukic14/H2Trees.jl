module H2BEASTTrees
using StaticArrays
using BEAST
using H2Trees
using Graphs

import H2Trees: TwoNTree, boundingbox, numberoflevels, isgalerkinsymmetric
import H2Trees: BEASTProtrusionFunctor, halfsize, center, adjacencygraph
import H2Trees: MetisForest, MetisTree, noboundingsphereupdate
"""
    TwoNTree(space::BEAST.Space, minhalfsize; kwargs...)

Construct a TwoNTree from a given `BEAST.Space`.

# Arguments

  - `space::BEAST.Space`: The input space.
  - `minhalfsize`: The minimum half-size of the tree.
  - `computeprotrusion`: Protrusion functor used during tree construction. Defaults to `BEASTProtrusionFunctor(space)`.
  - `kwargs...`: Additional keyword arguments.

# Returns

A TwoNTree.
"""
function TwoNTree(
    space::BEAST.Space,
    minhalfsize;
    computeprotrusion=BEASTProtrusionFunctor(space),
    kwargs...,
)
    return TwoNTree(
        BEAST.positions(space), minhalfsize; computeprotrusion=computeprotrusion, kwargs...
    )
end

"""
    TwoNTree(testspace::BEAST.Space, trialspace::BEAST.Space, minhalfsize; kwargs...)

Construct a block tree with two `TwoNTree`s from two given spaces: a test space and a
trial space.

# Arguments

  - `testspace::BEAST.Space`: The test space.
  - `trialspace::BEAST.Space`: The trial space.
  - `minhalfsize`: The minimum half-size of the tree.
  - `testcomputeprotrusion`: Protrusion functor for the test tree. Defaults to `BEASTProtrusionFunctor(testspace)`.
  - `trialcomputeprotrusion`: Protrusion functor for the trial tree. Defaults to `BEASTProtrusionFunctor(trialspace)`.
  - `kwargs...`: Additional keyword arguments.

# Returns

A TwoNTree.
"""
function TwoNTree(
    testspace::BEAST.Space,
    trialspace::BEAST.Space,
    minhalfsize;
    testcomputeprotrusion=BEASTProtrusionFunctor(testspace),
    trialcomputeprotrusion=BEASTProtrusionFunctor(trialspace),
    kwargs...,
)
    return TwoNTree(
        BEAST.positions(testspace),
        BEAST.positions(trialspace),
        minhalfsize;
        testcomputeprotrusion=testcomputeprotrusion,
        trialcomputeprotrusion=trialcomputeprotrusion,
        kwargs...,
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
