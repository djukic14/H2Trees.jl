using Test
using ParallelKMeans
using PlotlyJS
using H2Trees
using StaticArrays
using CompScienceMeshes

@testset "H2PlotlyJSTrees" begin
    @test_nowarn traceball(SVector(0.0, 0.0, 0.0), 1.0; n=30)
    @test typeof(traceball(SVector(0.0, 0.0, 0.0), 1.0; n=30)) ==
        GenericTrace{Dict{Symbol,Any}}
    @test sort(collect(keys(traceball(SVector(0.0, 0.0, 0.0), 1.0; n=30).fields))) ==
        [:type, :x, :y, :z]

    @test_nowarn traceball(SVector(0.0, 0.0), 1.0; n=30)
    @test typeof(traceball(SVector(0.0, 0.0), 1.0; n=30)) == GenericTrace{Dict{Symbol,Any}}
    @test sort(collect(keys(traceball(SVector(0.0, 0.0), 1.0; n=30).fields))) ==
        [:type, :x, :y]
    @test_nowarn tracecube(SVector(0.0, 0.0, 0.0), 1.0)
    @test typeof(tracecube(SVector(0.0, 0.0, 0.0), 1.0)) == GenericTrace{Dict{Symbol,Any}}
    @test sort(collect(keys(tracecube(SVector(0.0, 0.0, 0.0), 1.0).fields))) ==
        [:type, :x, :y, :z]

    @test_nowarn tracecube(SVector(0.0, 0.0), 1.0)
    @test typeof(tracecube(SVector(0.0, 0.0), 1.0)) == GenericTrace{Dict{Symbol,Any}}
    @test sort(collect(keys(tracecube(SVector(0.0, 0.0), 1.0).fields))) == [:type, :x, :y]

    m = CompScienceMeshes.readmesh(
        joinpath(pkgdir(H2Trees), "test", "assets", "in", "multiplerects2.in")
    )
    tree = TwoNTree(vertices(m), 0.5)
    @test_nowarn tracecube(tree, 4)
    @test typeof(tracecube(tree, 4)) == GenericTrace{Dict{Symbol,Any}}
    @test sort(collect(keys(tracecube(tree, 4)))) == [:type, :x, :y, :z]

    tree = KMeansTree(vertices(m), 10; minvalues=100, n_threads=1)
    @test_nowarn traceball(tree, 4)
    @test typeof(traceball(tree, 4)) == GenericTrace{Dict{Symbol,Any}}
    @test sort(collect(keys(traceball(tree, 4)))) == [:type, :x, :y, :z]
end
