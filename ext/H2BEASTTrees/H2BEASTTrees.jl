module H2BEASTTrees
using StaticArrays
using BEAST
using H2Trees

import H2Trees: TwoNTree, QuadPointsTree, boundingbox, numberoflevels, isgalerkinsymmetric

"""
    TwoNTree(space::BEAST.Space, minhalfsize; kwargs...)

Construct a TwoNTree from a given `BEAST.Space`.

# Arguments

  - `space::BEAST.Space`: The input space.
  - `minhalfsize`: The minimum half-size of the tree.
  - `kwargs...`: Additional keyword arguments.

# Returns

A TwoNTree.
"""
function TwoNTree(space::BEAST.Space, minhalfsize; kwargs...)
    return TwoNTree(BEAST.positions(space), minhalfsize; kwargs...)
end

"""
    TwoNTree(testspace::BEAST.Space, trialspace::BEAST.Space, minhalfsize; kwargs...)

Construct a block tree with two `TwoNTree`s from two given spaces: a test space and a
trial space.

# Arguments

  - `testspace::BEAST.Space`: The test space.
  - `trialspace::BEAST.Space`: The trial space.
  - `minhalfsize`: The minimum half-size of the tree.
  - `kwargs...`: Additional keyword arguments.

# Returns

A TwoNTree.
"""
function TwoNTree(testspace::BEAST.Space, trialspace::BEAST.Space, minhalfsize; kwargs...)
    return TwoNTree(
        BEAST.positions(testspace), BEAST.positions(trialspace), minhalfsize; kwargs...
    )
end

include("QuadPointsTree.jl")

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

end # module H2BEASTTrees
