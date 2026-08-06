using CompScienceMeshes
using H2Trees
using PlotlyJS
using ParallelKMeans

m = meshsphere(1.0, 0.1)
tree = KMeansTree(
    vertices(m); builder=KMeansTreeBuilder(; numberofclusters=4, minvalues=60)
)

traces = [wireframe(skeleton(m, 1))]

for node in H2Trees.LevelIterator(tree, 4)
    push!(
        traces,
        H2Trees.traceball(
            tree, node; colorscale=[[0, :pink], [1, :pink]], opacity=0.6, showscale=false
        ),
    )
end

p = PlotlyJS.plot(
    traces,
    Layout(;
        scene=attr(;
            xaxis=attr(; visible=false),
            yaxis=attr(; visible=false),
            zaxis=attr(; visible=false),
        ),
    ),
)
