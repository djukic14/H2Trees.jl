using Test
using PlotlyJS
using H2Trees
using StaticArrays

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
end
