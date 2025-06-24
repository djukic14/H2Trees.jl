module H2BEASTTrees
using StaticArrays
using BEAST
using H2Trees

import H2Trees: TwoNTree, QuadPointsTree, boundingbox, numberoflevels, isgalerkinsymmetric

function TwoNTree(space::BEAST.Space, minhalfsize; kwargs...)
    return TwoNTree(BEAST.positions(space), minhalfsize; kwargs...)
end

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
