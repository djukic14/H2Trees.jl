using Test
using Metis, Graphs, LinearAlgebra, SparseArrays
using CompScienceMeshes, BEAST
using H2Trees

function dual_piecewise_constant(
    mesh,
    jct=skeleton(m, 0),# jct=CompScienceMeshes.mesh(coordtype(mesh), dimension(mesh) - 1)
)
    vertexlist = BEAST.interior_and_junction_vertices(mesh, jct)
    return dual_piecewise_constant(mesh, vertexlist)
end

function dual_piecewise_constant(mesh, vertexlist::Vector;)
    T = coordtype(mesh)

    fns = Vector{Vector{BEAST.Shape{T}}}(undef, length(vertexlist))
    pos = Vector{vertextype(mesh)}()

    fine = barycentric_refinement(mesh)
    vtoc, vton = vertextocellmap(fine)
    for (k, v) in enumerate(vertexlist)
        n = vton[v]
        F = vtoc[v, 1:n]
        fns[k] = single_dual_piecewise_constant(fine, F, v)
        push!(pos, mesh.vertices[v])
    end

    NF = 1
    return BEAST.LagrangeBasis{0,-1,NF}(fine, fns, pos)
end

function single_dual_piecewise_constant(fine, F, v; normalize=false)
    T = coordtype(fine)
    fn = BEAST.Shape{T}[]

    area = T(0)

    for cellid in F
        ptch = chart(fine, cellid)
        area += CompScienceMeshes.volume(ptch)
    end

    for cellid in F
        coeff = 1.0 / area
        refid = 1
        push!(fn, BEAST.Shape(cellid, refid, coeff))
    end

    return fn
end

m = meshicosphere(1)

# m = meshrectangle(1.0, 1.0, 1.0)
# m = barycentric_refinement(m).mesh

verts = CompScienceMeshes.skeleton(m, 0)
edges = CompScienceMeshes.skeleton(m, 1)
Λ = CompScienceMeshes.connectivity(verts, edges)

# edges = setminus(skeleton(m, 1), boundary(m))
# verts = setminus(skeleton(m, 0), skeleton(boundary(m), 0))
# Λ2 = CompScienceMeshes.connectivity(verts, edges, sign)

ΛΛ = Λ' * Λ
A = -ΛΛ
A[diagind(A)] .= 0

X = dual_piecewise_constant(m)
# X = lagrangec0d1(m)
Y = lagrangec0d1(m)

edges = BEAST.skeleton(geometry(X), 1)
faces = BEAST.skeleton(geometry(X), 2)
Σ = BEAST.connectivity(faces, edges, sign)
ΣΣ = Σ' * Σ

g, w = H2Trees.adjacencygraph(X);
gA = Graph(A)

tree = H2Trees.MetisTree(X, 4);

areas = [CompScienceMeshes.volume(chart(X.geo, i)) for i in 1:numcells(X.geo)]
functionarea = zeros(length(areas))

nodesarea = zeros(H2Trees.numberofnodes(tree))
for node in H2Trees.DepthFirstIterator(tree)
    nodesarea[node] = sum(w[H2Trees.values(tree, node)])
end

levelareas = [
    (nodesarea[H2Trees.LevelIterator(tree, level)]) for level in H2Trees.levels(tree)
]
maxminlevelareas = [
    (maximum(levelareas[i]), minimum(levelareas[i])) for i in eachindex(levelareas)
]

# sort!.(g.fadjlist)
# sort!.(gA.fadjlist)
#
# sort!(g.fadjlist)
# sort!(gA.fadjlist)

# using Plots, GraphRecipes

# graphplot(
#     g;
#     curves=false,
#     nodelabel=string.(1:nv(g)),
#     nodelabelsize=20,     # bigger
#     nodelabeldist=2.0,    # further away
#     nodesize=0.2,          # shrink nodes so labels aren't covered
# )

# using PlotlyJS

# function plotmeshlabels(m; labels=string.(1:numcells(m)))
#     centers = [cartesian(CompScienceMeshes.center(chart(m, i))) for i in 1:numcells(m)]
#     cx = [c[1] for c in centers]
#     cy = [c[2] for c in centers]
#     cz = [c[3] for c in centers]

#     label_trace = scatter3d(;
#         x=cx,
#         y=cy,
#         z=cz,
#         mode="text",
#         text=labels,
#         textposition="middle center",
#         textfont=attr(; color="red", size=10),
#     )

#     return [wireframe(skeleton(m, 1)), label_trace]
# end

# labels = Any["" for _ in 1:numcells(geometry(X))]

# for fns in sparse(g)[1, :].nzind
#     for fn in X.fns[fns]
#         labels[fn.cellid] = string(fn.cellid)
#     end
# end

# plot(plotmeshlabels(X.geo; labels=labels))

# # sparse(g)
# #
