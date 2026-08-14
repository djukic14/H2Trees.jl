# --- hide-from-docs ---
# (already shown verbatim in the first example block on this page)
using CompScienceMeshes, BEAST
using H2Trees

mx = meshsphere(1.0, 0.1)
my = meshsphere(2.0, 0.1)
translate!(my, [6.0, 0.0, 0.0])

tree = buildtree(
    lagrangecxd0(mx),
    lagrangecxd0(my);
    builder=BlockTreeBuilder(;
        test=TwoNTreeBuilder(; minhalfsize=0.1), trial=TwoNTreeBuilder(; minhalfsize=0.1)
    ),
)
# --- end-hide-from-docs ---

using PlotlyJS
traces = [wireframe(skeleton(mx, 1)), wireframe(skeleton(my, 1))]
for node in H2Trees.LevelIterator(H2Trees.testtree(tree), 4)
    push!(
        traces,
        H2Trees.tracecube(H2Trees.testtree(tree), node; mode="lines", line_color=:pink),
    )
end
for node in H2Trees.LevelIterator(H2Trees.trialtree(tree), 4)
    push!(
        traces,
        H2Trees.tracecube(H2Trees.trialtree(tree), node; mode="lines", line_color=:blue),
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
        showlegend=false,
    ),
)
