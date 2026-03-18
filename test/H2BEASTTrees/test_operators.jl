using Test
using BEAST, H2Trees

@testset "BEAST symmetric operators" begin
    k = 1.0
    op = Helmholtz3D.singlelayer(; wavenumber=k)
    @test H2Trees.isgalerkinsymmetric(op)

    op = Helmholtz3D.doublelayer(; wavenumber=k)
    @test !H2Trees.isgalerkinsymmetric(op)

    op = Helmholtz3D.hypersingular(; wavenumber=k)
    @test H2Trees.isgalerkinsymmetric(op)

    op = Maxwell3D.singlelayer(; wavenumber=k)
    @test H2Trees.isgalerkinsymmetric(op)

    op = Maxwell3D.doublelayer(; wavenumber=k)
    @test H2Trees.isgalerkinsymmetric(op)
end
