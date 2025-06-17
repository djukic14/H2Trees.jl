using H2Trees
using Test
using Aqua
using JET

@testset "H2Trees.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(H2Trees)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(H2Trees; target_defined_modules = true)
    end
    # Write your tests here.
end
