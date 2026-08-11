# Phase 1: Ball construction, validation, and public API surface.

@testset "Ball construction and validation" begin
    @testset "dimensions 1, 2, 3 accepted" begin
        @test SB.Ball(SVector(0.0), 1.0) isa SB.Ball{1,Float64}
        @test SB.Ball(SVector(0.0, 0.0), 1.0) isa SB.Ball{2,Float64}
        @test SB.Ball(SVector(0.0, 0.0, 0.0), 1.0) isa SB.Ball{3,Float64}
    end

    @testset "dimension 0 and 4 rejected" begin
        @test_throws ArgumentError SB.Ball(SVector{0,Float64}(), 1.0)
        @test_throws ArgumentError SB.Ball(SVector(0.0, 0.0, 0.0, 0.0), 1.0)
    end

    @testset "negative radius rejected, zero accepted" begin
        @test_throws ArgumentError SB.Ball(SVector(0.0, 0.0), -1.0)
        @test SB.radius(SB.Ball(SVector(0.0, 0.0), 0.0)) == 0.0
    end

    @testset "non-finite input rejected" begin
        @test_throws ArgumentError SB.Ball(SVector(NaN, 0.0), 1.0)
        @test_throws ArgumentError SB.Ball(SVector(Inf, 0.0), 1.0)
        @test_throws ArgumentError SB.Ball(SVector(0.0, 0.0), NaN)
        @test_throws ArgumentError SB.Ball(SVector(0.0, 0.0), Inf)
    end

    @testset "integer input promotion" begin
        b = SB.Ball([1, 2, 3], 4)
        @test b isa SB.Ball{3,Float64}
        @test SB.center(b) == SVector(1.0, 2.0, 3.0)
        @test SB.radius(b) == 4.0
    end

    @testset "mixed scalar promotion" begin
        b = SB.Ball(SVector(1.0f0, 2.0f0), 1)
        @test b isa SB.Ball{2,Float32}
    end

    @testset "accessors and eltype" begin
        b = SB.Ball(SVector(1.0, 2.0), 3.0)
        @test SB.center(b) === SVector(1.0, 2.0)
        @test SB.radius(b) === 3.0
        @test eltype(typeof(b)) === Float64
        @test SB.dimension(b) == 2
    end
end

@testset "smallest_enclosing_ball input validation" begin
    @testset "empty input rejected" begin
        @test_throws ArgumentError SB.smallest_enclosing_ball(SB.Ball{2,Float64}[])
    end

    @testset "centers/radii length mismatch" begin
        @test_throws DimensionMismatch SB.smallest_enclosing_ball(
            [SVector(0.0, 0.0)], [1.0, 2.0]
        )
    end

    @testset "output type preservation" begin
        b64 = [SB.Ball(SVector(0.0, 0.0), 1.0), SB.Ball(SVector(2.0, 0.0), 1.0)]
        @test SB.smallest_enclosing_ball(b64) isa SB.Ball{2,Float64}
        b32 = [SB.Ball(SVector(0.0f0, 0.0f0), 1.0f0), SB.Ball(SVector(2.0f0, 0.0f0), 1.0f0)]
        @test SB.smallest_enclosing_ball(b32) isa SB.Ball{2,Float32}
    end

    @testset "centers/radii convenience method" begin
        centers = [SVector(0.0, 0.0), SVector(2.0, 0.0)]
        radii = [1.0, 1.0]
        r = SB.smallest_enclosing_ball(centers, radii)
        @test r isa SB.Ball{2,Float64}
        @test balls_approx(r, SB.Ball(SVector(1.0, 0.0), 2.0))
    end
end

@testset "encloses predicate (independent residual cross-check)" begin
    outer = SB.Ball(SVector(0.0, 0.0), 2.0)
    inner_in = SB.Ball(SVector(0.5, 0.0), 1.0)
    inner_out = SB.Ball(SVector(3.0, 0.0), 1.0)
    @test SB.encloses(outer, inner_in)
    @test !SB.encloses(outer, inner_out)
    # Cross-check the boolean against an independently computed residual.
    @test (indep_residual(outer, inner_in) <= 1e-12) == SB.encloses(outer, inner_in)
    @test (indep_residual(outer, inner_out) <= 1e-12) == SB.encloses(outer, inner_out)
end
