using Test

@testset verbose = true "1D" begin
    include("test_1d.jl")
end

@testset verbose = true "2D" begin
    include("test_2d.jl")
end

@testset verbose = true "3D" begin
    include("test_3d.jl")
end

@testset verbose = true "Properties" begin
    include("test_properties.jl")
end

@testset verbose = true "Tree layout" begin
    include("test_treelayout.jl")
end
