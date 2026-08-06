using H2Trees
using Metis
using BEAST
using CompScienceMeshes
using PlotlyJS

m = meshicosphere(10)
X = lagrangecxd0(m)

tree = H2Trees.MetisTree(X; builder=H2Trees.MetisTreeBuilder(; numdivisions=4))

childcolors = [
    "rgb(0.5490196078431373, 0.00784313725490196, 0.45098039215686275)",
    "rgb(0.6, 0.38823529411764707, 0.18823529411764706)",
    "rgb(0.5019607843137255, 0.7725490196078432, 0.37254901960784315)",
    "rgb(0.7019607843137254, 0.9490196078431372, 0.9921568627450981)",
]

color = ["rgba(0,0,0,0)" for _ in 1:numcells(m)]

node = H2Trees.root(tree)
for (i, child) in enumerate(H2Trees.children(tree, node))
    for val in H2Trees.values(tree, child)
        for fn in X.fns[val]
            color[fn.cellid] = childcolors[i]
        end
    end
end

# --- hide-from-docs ---
traces = [wireframe(skeleton(m, 1); width=1)]
push!(traces, mesh3d(m; facecolor=color))

p = PlotlyJS.plot(
    traces,
    Layout(;
        scene=attr(;
            aspectmode="data",
            xaxis=attr(; visible=false),
            yaxis=attr(; visible=false),
            zaxis=attr(; visible=false),
        ),
    ),
)
