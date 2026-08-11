# Phase 5: degeneracy, near-degeneracy, and precision robustness.

@testset "duplicate balls" begin
    b = SB.Ball(SVector(1.0, 2.0), 1.5)
    r = SB.smallest_enclosing_ball([b, b, b])
    @test balls_approx(r, b)
end

@testset "duplicate centers, unequal radii" begin
    balls = [
        SB.Ball(SVector(0.0, 0.0), 1.0),
        SB.Ball(SVector(0.0, 0.0), 2.5),
        SB.Ball(SVector(0.0, 0.0), 0.3),
    ]
    r = SB.smallest_enclosing_ball(balls)
    @test balls_approx(r, SB.Ball(SVector(0.0, 0.0), 2.5))
end

@testset "several identical largest balls" begin
    big = SB.Ball(SVector(0.0, 0.0, 0.0), 4.0)
    balls = [big, big, big, SB.Ball(SVector(1.0, 0.0, 0.0), 0.5)]
    r = SB.smallest_enclosing_ball(balls)
    @test balls_approx(r, big)
end

@testset "collinear centers (2D and 3D)" begin
    b2 = [SB.Ball(SVector(Float64(i), 0.0), 0.5) for i in 0:4]
    r2 = SB.smallest_enclosing_ball(b2)
    @test balls_approx(r2, SB.Ball(SVector(2.0, 0.0), 2.5); rtol=1e-7)
    assert_encloses(r2, b2)

    b3 = [SB.Ball(SVector(Float64(i), Float64(i), Float64(i)), 0.2) for i in 0:3]
    r3 = SB.smallest_enclosing_ball(b3)
    assert_encloses(r3, b3)
    # collinear 3D -> answer is a two-ball diameter
    cert = SB._smallest_enclosing_ball_with_certificate(b3)
    @test cert.support_size <= 2
end

@testset "coplanar centers in 3D" begin
    pts = [
        SVector(1.0, 0.0, 0.0),
        SVector(-1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, -1.0, 0.0),
    ]
    balls = [SB.Ball(p, 0.3) for p in pts]
    r = SB.smallest_enclosing_ball(balls)
    @test balls_approx(r, SB.Ball(SVector(0.0, 0.0, 0.0), 1.3); rtol=1e-7)
    assert_encloses(r, balls)
    cert = SB._smallest_enclosing_ball_with_certificate(balls)
    @test cert.support_size <= 3   # coplanar square: 2 antipodal balls suffice actually
end

@testset "near-collinear / near-coplanar do not crash and enclose" begin
    b = [
        SB.Ball(SVector(0.0, 0.0), 0.5),
        SB.Ball(SVector(1.0, 1e-9), 0.5),
        SB.Ball(SVector(2.0, 0.0), 0.5),
    ]
    r = SB.smallest_enclosing_ball(b)
    assert_encloses(r, b)
end

@testset "very small nonzero radii" begin
    balls = [SB.Ball(SVector(0.0, 0.0), 1e-12), SB.Ball(SVector(1e-11, 0.0), 1e-12)]
    r = SB.smallest_enclosing_ball(balls)
    assert_encloses(r, balls)
end

@testset "radii differing by many orders of magnitude" begin
    balls = [
        SB.Ball(SVector(0.0, 0.0), 1e8),
        SB.Ball(SVector(1.0, 0.0), 1e-6),
        SB.Ball(SVector(-2.0, 3.0), 1.0),
    ]
    r = SB.smallest_enclosing_ball(balls)
    @test balls_approx(r, balls[1]; rtol=1e-6)
    assert_encloses(r, balls)
end

@testset "large translation, small relative separation" begin
    base = SVector(1e6, 1e6)
    balls = [SB.Ball(base, 1.0), SB.Ball(base .+ SVector(3.0, 0.0), 1.0)]
    r = SB.smallest_enclosing_ball(balls)
    @test balls_approx(r, SB.Ball(base .+ SVector(1.5, 0.0), 2.5); rtol=1e-8)
end

@testset "exact vs near internal tangency transitions" begin
    # outer nominally R=5 centered origin, inner at distance 4 radius 1 -> exactly tangent.
    for ε in (-1e-12, -1e-14, 0.0, 1e-14, 1e-12)
        balls = [SB.Ball(SVector(0.0, 0.0), 5.0), SB.Ball(SVector(4.0 + ε, 0.0), 1.0)]
        r = SB.smallest_enclosing_ball(balls)
        assert_encloses(r, balls; rtol=1e-9)
        # radius should be ~5 when inner is inside/tangent, slightly more when poking out
        @test SB.radius(r) >= 5.0 - 1e-9
    end
end

@testset "symmetric instance with several equivalent supports" begin
    # regular hexagon vertices, equal radii -> unique center at origin
    pts = [SVector(cos(2π * k / 6), sin(2π * k / 6)) for k in 0:5]
    balls = [SB.Ball(p, 0.25) for p in pts]
    r = SB.smallest_enclosing_ball(balls)
    @test balls_approx(r, SB.Ball(SVector(0.0, 0.0), 1.25); rtol=1e-7)
    # permutation invariance under this symmetry
    rng = MersenneTwister(99)
    for _ in 1:10
        rp = SB.smallest_enclosing_ball(balls[randperm(rng, 6)])
        @test balls_approx(rp, r)
    end
end

@testset "mesh-derived near-collinear triple + centered ball (regression)" begin
    # Real MetisTree child balls (spherewithcenter mesh) that made the SEBB enumeration crash
    # with "failed to find a numerically valid enclosing ball". Three zero-radius point balls
    # are *exactly* collinear (b1,b2,b3 with b2 the midpoint), and a fourth ball of radius
    # ~0.049 sits a hair off that line near the midpoint. The optimal support is the
    # near-affinely-dependent triple {b1, b3, b4}, whose tangent-ball construction hits a
    # near-double quadratic root, so the tangency/containment residuals are limited by the
    # support's conditioning (~1e-7) rather than machine precision. Every candidate used to be
    # rejected; the conditioning-aware tangency/containment slack now accepts it.
    balls = [
        SB.Ball(SVector(-0.594389364929929, 0.35174161864365217, -0.7203936662983679), 0.0),
        SB.Ball(
            SVector(-0.5698313417744849, 0.39380908116689756, -0.7185003068964156), 0.0
        ),
        SB.Ball(
            SVector(-0.5439886765789463, 0.43497390667837965, -0.7147829492191635), 0.0
        ),
        SB.Ball(
            SVector(-0.5693832694025134, 0.39339011800138146, -0.7174609659287856),
            0.04864210473407279,
        ),
    ]
    r = SB.smallest_enclosing_ball(balls)
    @test r isa SB.Ball{3,Float64}
    assert_encloses(r, balls; rtol=1e-6)
    # The result must be no larger than the enclosing radius of the biggest child alone would
    # need, and at least the largest child's radius.
    @test SB.radius(r) >= 0.04864210473407279 - 1e-9
    @test SB.radius(r) <= 0.05

    # Permutation invariance despite the degeneracy.
    rng = MersenneTwister(2024)
    for _ in 1:16
        rp = SB.smallest_enclosing_ball(balls[randperm(rng, 4)])
        @test balls_approx(rp, r; rtol=1e-6)
    end

    # Certificate: convex weights over an active support of size <= 4.
    cert = SB._smallest_enclosing_ball_with_certificate(balls)
    @test 1 <= cert.support_size <= 4
    w = SB.active_weights(cert)
    @test all(>=(-1e-9), w)
    @test isapprox(sum(w), 1.0; atol=1e-7)
end

@testset "Float32" begin
    balls = [
        SB.Ball(SVector(0.0f0, 0.0f0), 1.0f0),
        SB.Ball(SVector(3.0f0, 0.0f0), 0.5f0),
        SB.Ball(SVector(1.0f0, 2.0f0), 0.8f0),
    ]
    r = SB.smallest_enclosing_ball(balls)
    @test r isa SB.Ball{2,Float32}
    assert_encloses(r, balls; rtol=1.0f-4)
end

@testset "BigFloat matches Float64" begin
    setprecision(BigFloat, 256) do
        pts = [SVector(cos(2π * k / 3), sin(2π * k / 3)) for k in 0:2]
        balls64 = [SB.Ball(p, 0.4) for p in pts]
        ptsbig = [
            SVector(BigFloat(cos(2 * big(π) * k / 3)), BigFloat(sin(2 * big(π) * k / 3)))
            for k in 0:2
        ]
        ballsbig = [SB.Ball(p, big"0.4") for p in ptsbig]
        r64 = SB.smallest_enclosing_ball(balls64)
        rbig = SB.smallest_enclosing_ball(ballsbig)
        @test rbig isa SB.Ball{2,BigFloat}
        @test abs(Float64(SB.radius(rbig)) - SB.radius(r64)) <= 1e-7
        @test norm(Float64.(SB.center(rbig)) .- SB.center(r64)) <= 1e-7
    end
end
