using H2Trees
using Metis
using BEAST
using CompScienceMeshes
using PlotlyJS

m = CompScienceMeshes.readmesh(
    joinpath(pkgdir(H2Trees), "test", "assets", "in", "twospheres2.in")
)
X = lagrangecxd0(m)

forest = H2Trees.MetisForest(
    X;
    builder=H2Trees.MetisForestBuilder(;
        treebuilder=H2Trees.MetisTreeBuilder(; numdivisions=4)
    ),
)

childcolors = [
    "rgb(0.5490196078431373, 0.00784313725490196, 0.45098039215686275)",
    "rgb(0.6, 0.38823529411764707, 0.18823529411764706)",
    "rgb(0.5019607843137255, 0.7725490196078432, 0.37254901960784315)",
    "rgb(0.7019607843137254, 0.9490196078431372, 0.9921568627450981)",
    "rgb(0.5725490196078431, 0.16470588235294117, 0.34901960784313724)",
    "rgb(0.615686274509804, 0.5137254901960784, 0.11764705882352941)",
    "rgb(0.4, 0.8470588235294118, 0.611764705882353)",
    "rgb(0.4235294117647059, 0.9215686274509803, 0.8588235294117647)",
]

color = ["rgba(0,0,0,0)" for _ in 1:numcells(m)]

node = H2Trees.root(forest[1])
for (i, child) in enumerate(H2Trees.children(forest[1], node))
    for val in H2Trees.values(forest[1], child)
        for fn in X.fns[val]
            color[fn.cellid] = childcolors[i]
        end
    end
end

node = H2Trees.root(forest[2])
for (i, child) in enumerate(H2Trees.children(forest[2], node))
    for val in H2Trees.values(forest[2], child)
        for fn in X.fns[val]
            color[fn.cellid] = childcolors[i + 4]
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
