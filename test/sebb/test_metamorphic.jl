# Phase 4: metamorphic tests (transformations that predict how the answer must change).

# A fixed, moderately generic set of base instances used across transforms.
function base_instances_2d()
    return [
        [
            SB.Ball(SVector(0.0, 0.0), 1.0),
            SB.Ball(SVector(3.0, 0.0), 0.5),
            SB.Ball(SVector(1.0, 2.0), 0.8),
        ],
        [SB.Ball(SVector(cos(2π * k / 3), sin(2π * k / 3)), 0.4) for k in 0:2],
        [
            SB.Ball(SVector(-1.0, -1.0), 0.2),
            SB.Ball(SVector(2.0, 1.0), 1.3),
            SB.Ball(SVector(0.0, 3.0), 0.7),
            SB.Ball(SVector(0.5, 0.5), 0.1),
        ],
    ]
end

function base_instances_3d()
    return [
        [
            SB.Ball(SVector(1.0, 1.0, 1.0), 0.3),
            SB.Ball(SVector(1.0, -1.0, -1.0), 0.3),
            SB.Ball(SVector(-1.0, 1.0, -1.0), 0.3),
            SB.Ball(SVector(-1.0, -1.0, 1.0), 0.3),
        ],
        [
            SB.Ball(SVector(0.0, 0.0, 0.0), 1.0),
            SB.Ball(SVector(2.0, 0.0, 0.0), 0.5),
            SB.Ball(SVector(0.0, 2.0, 1.0), 0.8),
            SB.Ball(SVector(1.0, 1.0, 3.0), 0.2),
        ],
    ]
end

@testset "permutation invariance" begin
    rng = MersenneTwister(1)
    for balls in vcat(base_instances_2d(), base_instances_3d())
        ref = SB.smallest_enclosing_ball(balls)
        for _ in 1:20
            p = randperm(rng, length(balls))
            r = SB.smallest_enclosing_ball(balls[p])
            @test balls_approx(r, ref)
        end
    end
end

@testset "translation equivariance" begin
    for balls in base_instances_2d()
        ref = SB.smallest_enclosing_ball(balls)
        for t in (SVector(10.0, -5.0), SVector(1e4, 1e4), SVector(-3.3, 2.1))
            tb = [SB.Ball(SB.center(b) .+ t, SB.radius(b)) for b in balls]
            r = SB.smallest_enclosing_ball(tb)
            @test isapprox(SB.radius(r), SB.radius(ref); rtol=1e-6)
            @test norm(SB.center(r) .- (SB.center(ref) .+ t)) <= 1e-5 * max(1.0, norm(t))
        end
    end
end

@testset "rotation/reflection equivariance" begin
    rng = MersenneTwister(2)
    # 2D rotations
    for balls in base_instances_2d()
        ref = SB.smallest_enclosing_ball(balls)
        for θ in (0.3, 1.7, -2.2)
            Q = SMatrix{2,2}(cos(θ), sin(θ), -sin(θ), cos(θ))
            rb = [SB.Ball(Q * SB.center(b), SB.radius(b)) for b in balls]
            r = SB.smallest_enclosing_ball(rb)
            @test isapprox(SB.radius(r), SB.radius(ref); rtol=1e-6)
            @test norm(SB.center(r) .- Q * SB.center(ref)) <=
                1e-6 * max(1.0, SB.radius(ref))
        end
        # reflection
        F = SMatrix{2,2}(1.0, 0.0, 0.0, -1.0)
        rb = [SB.Ball(F * SB.center(b), SB.radius(b)) for b in balls]
        r = SB.smallest_enclosing_ball(rb)
        @test norm(SB.center(r) .- F * SB.center(ref)) <= 1e-6 * max(1.0, SB.radius(ref))
    end
    # random orthogonal 3D matrices (via QR of random matrix)
    for balls in base_instances_3d()
        ref = SB.smallest_enclosing_ball(balls)
        for _ in 1:5
            A = SMatrix{3,3,Float64}(randn(rng, 3, 3))
            Q = Matrix(qr(A).Q)
            Qs = SMatrix{3,3}(Q)
            rb = [SB.Ball(Qs * SB.center(b), SB.radius(b)) for b in balls]
            r = SB.smallest_enclosing_ball(rb)
            @test isapprox(SB.radius(r), SB.radius(ref); rtol=1e-6)
            @test norm(SB.center(r) .- Qs * SB.center(ref)) <=
                1e-6 * max(1.0, SB.radius(ref))
        end
    end
end

@testset "positive scaling equivariance" begin
    for balls in vcat(base_instances_2d(), base_instances_3d())
        ref = SB.smallest_enclosing_ball(balls)
        for α in (0.001, 0.5, 2.0, 1e3)
            sb = [SB.Ball(α .* SB.center(b), α * SB.radius(b)) for b in balls]
            r = SB.smallest_enclosing_ball(sb)
            @test balls_approx(
                r, SB.Ball(α .* SB.center(ref), α * SB.radius(ref)); rtol=1e-6
            )
        end
    end
end

@testset "common radius shift" begin
    for balls in base_instances_2d()
        ref = SB.smallest_enclosing_ball(balls)
        for δ in (0.5, 2.0)
            sb = [SB.Ball(SB.center(b), SB.radius(b) + δ) for b in balls]
            r = SB.smallest_enclosing_ball(sb)
            @test norm(SB.center(r) .- SB.center(ref)) <= 1e-6 * max(1.0, SB.radius(ref))
            @test isapprox(SB.radius(r), SB.radius(ref) + δ; rtol=1e-6)
        end
    end
end

@testset "adding redundant (strictly contained) balls" begin
    for balls in base_instances_2d()
        ref = SB.smallest_enclosing_ball(balls)
        c = SB.center(ref)
        # a ball strictly inside ref
        extra = SB.Ball(c .+ SVector(0.01, 0.0), max(SB.radius(ref) - 0.5, 0.0) * 0.5)
        r = SB.smallest_enclosing_ball(vcat(balls, [extra]))
        @test balls_approx(r, ref)
    end
end

@testset "monotonicity: adding a ball cannot shrink the radius" begin
    rng = MersenneTwister(3)
    for balls in base_instances_3d()
        ref = SB.smallest_enclosing_ball(balls)
        for _ in 1:10
            extra = SB.Ball(SVector(randn(rng, 3)...) .* 3, abs(randn(rng)) * 0.5)
            r = SB.smallest_enclosing_ball(vcat(balls, [extra]))
            @test SB.radius(r) >= SB.radius(ref) - 1e-9
        end
    end
end

@testset "dimension embedding" begin
    # 1D -> 2D -> 3D by appending zeros must give the embedded result.
    c1 = [0.0, 3.0, -2.0, 5.0]
    r1 = [1.0, 0.5, 2.0, 0.3]
    balls1 = [SB.Ball(SVector(c1[i]), r1[i]) for i in eachindex(c1)]
    ref1 = SB.smallest_enclosing_ball(balls1)

    balls2 = [SB.Ball(SVector(c1[i], 0.0), r1[i]) for i in eachindex(c1)]
    ref2 = SB.smallest_enclosing_ball(balls2)
    @test isapprox(SB.radius(ref2), SB.radius(ref1); rtol=1e-7)
    @test isapprox(SB.center(ref2)[1], SB.center(ref1)[1]; atol=1e-7)
    @test abs(SB.center(ref2)[2]) <= 1e-7

    balls3 = [SB.Ball(SVector(c1[i], 0.0, 0.0), r1[i]) for i in eachindex(c1)]
    ref3 = SB.smallest_enclosing_ball(balls3)
    @test isapprox(SB.radius(ref3), SB.radius(ref1); rtol=1e-7)

    # 2D -> 3D (z=0 plane)
    for balls in base_instances_2d()
        ref = SB.smallest_enclosing_ball(balls)
        emb = [
            SB.Ball(SVector(SB.center(b)[1], SB.center(b)[2], 0.0), SB.radius(b)) for
            b in balls
        ]
        r = SB.smallest_enclosing_ball(emb)
        @test isapprox(SB.radius(r), SB.radius(ref); rtol=1e-6)
        @test abs(SB.center(r)[3]) <= 1e-6
    end
end
