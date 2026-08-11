# Phase 3/4: optimality certificates.
#
# Uses the qualified internal `_smallest_enclosing_ball_with_certificate`. For each result the
# certificate must be a convex combination (Lemma 2.2) of internally-tangent support balls
# (Lemma 3.1) that reproduces the center, with support size <= N + 1 (Lemma 2.4).

function check_certificate(balls; N)
    cert = SB._smallest_enclosing_ball_with_certificate(balls)
    ball = cert.ball
    T = typeof(SB.radius(ball))
    sup = SB.active_support(cert)
    w = SB.active_weights(cert)

    @test 1 <= cert.support_size <= N + 1
    @test length(sup) == cert.support_size
    @test allunique(sup)
    @test all(1 .<= sup .<= length(balls))

    # weights: nonnegative, sum to one
    wtol = sqrt(eps(T)) * 64
    @test all(>=(-wtol), w)
    @test abs(sum(w) - 1) <= wtol

    # weighted support centers reproduce the output center
    recon = sum(w[k] .* SB.center(balls[sup[k]]) for k in eachindex(sup))
    @test norm(recon .- SB.center(ball)) <= 1e-6 * max(1.0, norm(SB.center(ball)))

    # each support ball is internally tangent (independent residual ~ 0)
    for k in eachindex(sup)
        res = indep_residual(ball, balls[sup[k]])
        scale = max(SB.radius(ball), SB.radius(balls[sup[k]]), 1.0)
        @test abs(res) <= 1e-6 * scale
    end

    # every input enclosed
    assert_encloses(ball, balls)
    return cert
end

@testset "certificate: one-ball nested" begin
    big = SB.Ball(SVector(0.0, 0.0), 4.0)
    balls = [big, SB.Ball(SVector(1.0, 0.0), 0.5)]
    cert = check_certificate(balls; N=2)
    @test cert.support_size == 1
    @test SB.active_weights(cert) == [1.0]
end

@testset "certificate: two-ball" begin
    balls = [SB.Ball(SVector(0.0, 0.0), 1.0), SB.Ball(SVector(4.0, 0.0), 2.0)]
    check_certificate(balls; N=2)
end

@testset "certificate: triangle (three-ball support)" begin
    pts = [SVector(cos(2π * k / 3), sin(2π * k / 3)) for k in 0:2]
    balls = [SB.Ball(p, 0.4) for p in pts]
    cert = check_certificate(balls; N=2)
    @test cert.support_size == 3
end

@testset "certificate: tetrahedron (four-ball support)" begin
    tet = [
        SVector(1.0, 1.0, 1.0),
        SVector(1.0, -1.0, -1.0),
        SVector(-1.0, 1.0, -1.0),
        SVector(-1.0, -1.0, 1.0),
    ]
    balls = [SB.Ball(p, 0.2) for p in tet]
    cert = check_certificate(balls; N=3)
    @test cert.support_size == 4
end

@testset "certificate: support indices refer to ORIGINAL order" begin
    # Deliberately unsorted input; the biggest ball is at index 2. A nested result must cite 2.
    balls = [
        SB.Ball(SVector(1.0, 0.0), 0.3),
        SB.Ball(SVector(0.0, 0.0), 5.0),   # contains everything
        SB.Ball(SVector(-1.0, 1.0), 0.2),
    ]
    cert = SB._smallest_enclosing_ball_with_certificate(balls)
    @test cert.support_size == 1
    @test SB.active_support(cert) == [2]
end

@testset "certificate: analytic tangent-constructed instances" begin
    # Construct inputs tangent to a known output whose center is in the convex hull of
    # directions => certified optimal (Section 21).
    # 2D three-support, 120-degree directions.
    cstar = SVector(2.0, -1.0)
    Rstar = 3.0
    dirs = [SVector(cos(2π * k / 3), sin(2π * k / 3)) for k in 0:2]
    rs = [0.5, 1.0, 1.5]
    balls = [SB.Ball(cstar .- (Rstar - rs[i]) .* dirs[i], rs[i]) for i in 1:3]
    # add an interior ball strictly inside
    push!(balls, SB.Ball(cstar .+ SVector(0.1, 0.1), 0.2))
    cert = check_certificate(balls; N=2)
    @test balls_approx(cert.ball, SB.Ball(cstar, Rstar); rtol=1e-6)
end
