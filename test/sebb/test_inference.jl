# Phase 7: type-stability / inference.

@testset "inference: ball-vector method returns concrete Ball{N,T}" begin
    @testset "Float64" begin
        for N in (1, 2, 3)
            balls = [SB.Ball(SVector(ntuple(_ -> rand(), N)...), 0.5) for _ in 1:8]
            r = @inferred SB.smallest_enclosing_ball(balls)
            @test r isa SB.Ball{N,Float64}
        end
    end
    @testset "Float32" begin
        balls = [
            SB.Ball(SVector(rand(Float32), rand(Float32), rand(Float32)), 0.5f0) for
            _ in 1:5
        ]
        r = @inferred SB.smallest_enclosing_ball(balls)
        @test r isa SB.Ball{3,Float32}
    end
    @testset "varying counts" begin
        for n in (1, 2, 3, 4, 8)
            balls = [
                SB.Ball(SVector(Float64(i), Float64(i^2), Float64(-i)), 0.3) for i in 1:n
            ]
            r = @inferred SB.smallest_enclosing_ball(balls)
            @test r isa SB.Ball{3,Float64}
        end
    end
end

@testset "inference: certificate function" begin
    balls = [SB.Ball(SVector(Float64(i), 0.0, 0.0), 0.5) for i in 1:4]
    cert = @inferred SB._smallest_enclosing_ball_with_certificate(balls)
    @test cert isa SB.SEBBResult{3,Float64}
end

@testset "inference: predicates and roots" begin
    a = SB.Ball(SVector(0.0, 0.0), 2.0)
    b = SB.Ball(SVector(0.5, 0.0), 1.0)
    @test (@inferred SB.encloses(a, b)) isa Bool
    @test (@inferred SB.containment_residual(a, b)) isa Float64
    tol = SB.Tolerance{Float64}(0.0, 1e-8)
    @test (@inferred SB.real_roots(1.0, -3.0, 2.0, tol)) isa SB.QuadraticRoots{Float64}
end
