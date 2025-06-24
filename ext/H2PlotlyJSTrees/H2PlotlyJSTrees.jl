module H2PlotlyJSTrees
using PlotlyJS
using StaticArrays
using H2Trees
import H2Trees: center, radius, halfsize, H2ClusterTree, traceball, tracecube
include("tracecube.jl")
include("traceball.jl")
end # module H2PlotlyJSTrees
