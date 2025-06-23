module H2NURBSTrees

using NURBS
using BEAST
using H2Trees
using StaticArrays

import H2Trees: treetrait

const BEASTnurbs = Base.get_extension(BEAST, :BEASTnurbs)
if !isnothing(BEASTnurbs)
    include("QuadPointsTree.jl")
else
    @warn "Could not find BEASTnurbs extension"
end

end # H2NURBSTrees
