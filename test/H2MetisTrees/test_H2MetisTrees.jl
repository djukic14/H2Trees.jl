using Test
using Metis, Graphs, LinearAlgebra, SparseArrays
using CompScienceMeshes, BEAST
using H2Trees
H2MetisTrees = Base.get_extension(H2Trees, :H2MetisTrees)
H2BEASTTrees = Base.get_extension(H2Trees, :H2BEASTTrees)
meshes =
    [
        "cuboid",
        "cuboid2",
        "cuboid3",
        "cuboid4",
        "multiplerects",
        "multiplerects2",
        "multiplerects3",
        "sphere",
        "sphere2",
        "sphere3",
        "sphere4",
        "sphere5",
        "sphere6",
        "sphere7",
        "sphere8",
        "spherewithcenter",
        "spherewithcenter2",
        "spherewithcenter3",
        "spherewithcenter4",
        "spherewithcenter5",
        "spherewithcenter6",
        "spherewithcenter7",
        "spherewithcenter8",
        "spherewithcenter9",
        "spherewithcenter10",
        "spherewithcenter11",
        "spherewithcenter12",
        "spherewithcenter13",
        "spherewithcenter14",
        "spherewithcenter15",
        "spherewithcenter16",
        "twospheres",
        "twospheres2",
        "twospheres3",
        "twospheres4",
    ] .* ".in"

ms = [
    CompScienceMeshes.readmesh(joinpath(pkgdir(H2Trees), "test", "assets", "in", m)) for
    m in meshes
]

m = ms[5]
# m = meshicosphere(100)
# m = meshsphere(1.0, 0.1)
# m = meshrectangle(1.0, 1.0, 0.01)
X = lagrangecxd0(m)

edges = setminus(skeleton(X.geo, 1), boundary(X.geo))

areas = [volume(chart(X.geo, i)) for i in 1:numcells(X.geo)]
Σ = connectivity(X.geo, edges, sign);
ΣΣ = Σ' * Σ

# Σ2 = getstars(X.geo);

g, w = H2BEASTTrees.adjacencygraph(X);

A = -ΣΣ
A[diagind(A)] .= 0

# ΣΣ[diagindices(ΣΣ)] .= 0.0
gΣΣ = Graph(ΣΣ)
gA = Graph(A)

@test g == gA
@test maximum(abs, areas - w) < 1e-10
# g1 = Metis.graph(gΣΣ)
# g2 = Metis.graph(gA)

tree = H2MetisTrees.MetisTree(
    BEAST.positions(X), gΣΣ, areas, 4; minvalues=1, splitunconnectedpartitions=true
)

forest = H2Trees.MetisForest(X, 4; splitunconnectedpartitions=true)

# forest = H2Trees.MetisForest(raviartthomas(m), 4; splitunconnectedpartitions=true)

tree = H2Trees.MetisTree(X, 4; splitunconnectedpartitions=true)

nodeareas = [
    sum(areas[H2Trees.values(tree, node)]) for node in 1:H2Trees.numberofnodes(tree)
]

levelareas = [
    (nodeareas[H2Trees.LevelIterator(tree, level)]) for level in H2Trees.levels(tree)
]
maxminlevelareas = [
    (maximum(levelareas[i]), minimum(levelareas[i])) for i in eachindex(levelareas)
]
# gΣΣ = ΣΣ

# connectedparts = connected_components(gΣΣ)

# connectedparts[1] = [1, 2, 3, 4]

# gΣΣsub, localtoglobal = induced_subgraph(gΣΣ, connectedparts[1])

# parts = H2MetisTrees.metispartition(gΣΣsub, areas[localtoglobal], 4)

# partareas = zeros(length(unique(parts)))
# for i in eachindex(parts)
#     partareas[parts[i]] += areas[localtoglobal[i]]
# end

# partitions = [zeros(Int64, 0) for _ in 1:maximum(parts)]

# for (i, part) in enumerate(parts)
#     push!(partitions[part], i)
# end
# sort!.(partitions)

# @test length(unique(vcat(partitions...))) == size(ΣΣ, 1)

# using BoundingSphere

# center, radius = boundingsphere(X.pos)
# tree = H2Trees.BoundingBallTree
