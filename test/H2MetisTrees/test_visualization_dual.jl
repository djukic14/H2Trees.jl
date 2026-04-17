using Test
using Metis, Graphs, LinearAlgebra, SparseArrays
using CompScienceMeshes, BEAST
using H2Trees, PlotlyJS

childcolors = [
    "rgb(0.5490196078431373, 0.00784313725490196, 0.45098039215686275)",
    "rgb(0.5725490196078431, 0.16470588235294117, 0.34901960784313724)",
    "rgb(0.5882352941176471, 0.2784313725490196, 0.25882352941176473)",
    "rgb(0.6, 0.38823529411764707, 0.18823529411764706)",
    "rgb(0.615686274509804, 0.5137254901960784, 0.11764705882352941)",
    "rgb(0.592156862745098, 0.6627450980392157, 0.16470588235294117)",
    "rgb(0.5019607843137255, 0.7725490196078432, 0.37254901960784315)",
    "rgb(0.4, 0.8470588235294118, 0.611764705882353)",
    "rgb(0.4235294117647059, 0.9215686274509803, 0.8588235294117647)",
    "rgb(0.7019607843137254, 0.9490196078431372, 0.9921568627450981)",
][[1, 4, 7, 10, 2, 5, 8, 3, 6, 9]]

m = meshicosphere(10)
# X = dual_piecewise_constant(m)
X = duallagrangecxd0(m)
# X = lagrangec0d1(m)
# X = raviartthomas(m)

color = ["rgba(0,0,0,0)" for i in 1:numcells(geometry(X))]
traces = [wireframe(skeleton(m, 1); width=1)]
push!(traces, mesh3d(geometry(X); facecolor=color, opacity=0.9))

verts = CompScienceMeshes.skeleton(m, 0)
edges = CompScienceMeshes.skeleton(m, 1)
Λ = CompScienceMeshes.connectivity(verts, edges)
ΛΛ = Λ' * Λ

Σ = connectivity(m, edges, sign)

@show maximum(abs, Λ' * Σ)

g = Graph(ΛΛ)

tree = H2Trees.MetisTree(X, 4; minvalues=1)
# tree = H2Trees.MetisTree(
#     BEAST.positions(X),
#     g,
#     ones(numfunctions(X)),
#     4;
#     # minvalues=1,
#     # splitunconnectedpartitions=true,
# )

# node = 2
node = H2Trees.parent(tree, H2Trees.leaves(tree)[10]) # Get the parent of the first leaf
for (i, child) in enumerate(H2Trees.children(tree, node))
    for val in H2Trees.values(tree, child)
        for fn in X.fns[val]
            color[fn.cellid] = childcolors[i]
        end
        # color[val] = childcolors[i]
    end
end

plot(traces)
