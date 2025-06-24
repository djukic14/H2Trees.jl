module H2ParallelKMeansTrees
using ParallelKMeans
using H2Trees
import H2Trees: kmeanswrapper

function kmeanswrapper(points, numberofclusters::Int; kwargs...)
    return kmeans(points, numberofclusters; kwargs...)
end

end # module H2ParallelKMeansTrees
