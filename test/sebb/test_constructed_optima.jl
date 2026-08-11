# Constructed-optimum and support-size-transition tests.
#
# These tests build inputs whose expected answer is known WITHOUT running the SEBB
# implementation itself (per plan Section 21): a ball B(cstar, Rstar) tangent to every input
# ball is the true smallest enclosing ball exactly when the convex hull of the tangency
# directions contains the origin (Fischer & Gärtner 2004, Lemma 2.2). `tangent_balls` and
# `tetrahedron_dirs` (test/sebb/test_helpers.jl) build such instances.

@testset "3D constructed-optimum: tetrahedron support" begin
    cstar = SVector(2.0, -3.0, 0.5)
    Rstar = 7.0
    dirs = tetrahedron_dirs(Float64)
    radii = [0.2, 1.1, 2.0, 0.7]

    balls = tangent_balls(cstar, Rstar, dirs, radii)
    result = SB.smallest_enclosing_ball(balls)

    @test norm(SB.center(result) - cstar) <= 1e-10
    @test abs(SB.radius(result) - Rstar) <= 1e-10
    assert_encloses(result, balls)

    cert = SB._smallest_enclosing_ball_with_certificate(balls)
    @test cert.support_size == 4
    @test sort(SB.active_support(cert)) == [1, 2, 3, 4]
    w = SB.active_weights(cert)
    @test all(>=(-1e-10), w)
    @test isapprox(sum(w), 1.0; atol=1e-10)
end

@testset "3D constructed-optimum: with redundant interior balls" begin
    cstar = SVector(2.0, -3.0, 0.5)
    Rstar = 7.0
    dirs = tetrahedron_dirs(Float64)
    radii = [0.2, 1.1, 2.0, 0.7]
    support = tangent_balls(cstar, Rstar, dirs, radii)

    interior = [
        SB.Ball(cstar, 1.0),
        SB.Ball(cstar + SVector(0.3, -0.2, 0.1), 2.0),
        SB.Ball(cstar - SVector(0.5, 0.1, -0.2), 0.5),
    ]
    for b in interior
        @test norm(SB.center(b) - cstar) + SB.radius(b) < Rstar - 1e-8
    end
    balls = vcat(support, interior)

    result = SB.smallest_enclosing_ball(balls)
    @test norm(SB.center(result) - cstar) <= 1e-9
    @test abs(SB.radius(result) - Rstar) <= 1e-9
    assert_encloses(result, balls)
end

@testset "3D constructed-optimum: translated and scaled" begin
    cstar0 = SVector(2.0, -3.0, 0.5)
    Rstar0 = 7.0
    dirs = tetrahedron_dirs(Float64)
    radii0 = [0.2, 1.1, 2.0, 0.7]

    t = SVector(100.0, -50.0, 25.0)
    alpha = 3.5
    cstar = alpha * cstar0 + t
    Rstar = alpha * Rstar0
    radii = alpha .* radii0

    balls = tangent_balls(cstar, Rstar, dirs, radii)
    result = SB.smallest_enclosing_ball(balls)

    @test norm(SB.center(result) - cstar) <= 1e-7 * max(1.0, Rstar)
    @test abs(SB.radius(result) - Rstar) <= 1e-7 * max(1.0, Rstar)
    assert_encloses(result, balls)
end

@testset "3D constructed-optimum: unequal radii including one zero" begin
    cstar = SVector(-1.0, 4.0, 2.0)
    Rstar = 5.0
    dirs = tetrahedron_dirs(Float64)
    radii = [0.0, 0.5, 1.5, 3.0]

    balls = tangent_balls(cstar, Rstar, dirs, radii)
    result = SB.smallest_enclosing_ball(balls)

    @test norm(SB.center(result) - cstar) <= 1e-10
    @test abs(SB.radius(result) - Rstar) <= 1e-10
    assert_encloses(result, balls)
end

@testset "2D support-size transition: two-support base + moving third ball" begin
    # A two-ball support (opposite directions) already satisfies Lemma 2.2's convex-hull
    # criterion on its own (0.5 each on u and -u puts the origin in the hull), so it is the
    # true answer regardless of a third ball, unless that third ball is pushed far enough
    # out that it can no longer be enclosed without growing the ball.
    cstar = SVector(3.0, -2.0)
    Rstar = 5.0
    u = SVector(1.0, 0.0)
    v = SVector(0.0, 1.0)
    r1, r2, r3 = 1.0, 1.5, 0.8

    base = tangent_balls(cstar, Rstar, [u, -u], [r1, r2])

    for delta in (-0.5, 0.0, 0.5)
        # delta > 0: ball 3 sits strictly inside B(cstar, Rstar) (redundant).
        # delta < 0: ball 3 pokes strictly outside B(cstar, Rstar).
        # delta == 0: ball 3 is exactly tangent to B(cstar, Rstar).
        b3center = cstar - (Rstar - r3 - delta) * v
        balls = vcat(base, [SB.Ball(b3center, r3)])

        result = SB.smallest_enclosing_ball(balls)
        assert_encloses(result, balls)
        cert = SB._smallest_enclosing_ball_with_certificate(balls)

        if delta > 0
            @test balls_approx(result, SB.Ball(cstar, Rstar); rtol=1e-9)
            @test cert.support_size == 2
            @test 3 ∉ SB.active_support(cert)
        elseif delta < 0
            @test SB.radius(result) > Rstar + 1e-9
        else
            @test SB.radius(result) ≈ Rstar atol = 1e-8
            @test norm(SB.center(result) - cstar) <= 1e-8
        end
    end
end

@testset "3D support-size transition" begin
    cstar = SVector(0.0, 0.0, 0.0)
    Rstar = 5.0
    dirs = tetrahedron_dirs(Float64)
    radii = [0.5, 0.5, 0.5, 0.5]

    for delta in (-1e-8, 0.0, 1e-8)
        balls = tangent_balls(cstar, Rstar, dirs, radii)

        # Move fourth center along its tangent direction:
        # delta < 0 makes it interior, delta > 0 makes it slightly outside target.
        b4 = balls[4]
        balls[4] = SB.Ball(SB.center(b4) - delta * dirs[4], SB.radius(b4))

        result = SB.smallest_enclosing_ball(balls)
        assert_encloses(result, balls)

        cert = SB._smallest_enclosing_ball_with_certificate(balls)

        if delta < 0
            # |delta| here (1e-8) sits right at the default rtol ~ sqrt(eps(Float64)) ~ 1.5e-8
            # scale, so whether the fourth ball is resolved as strictly interior (support size
            # 3) or still-tangent-within-tolerance (support size 4) is legitimately
            # tolerance-dependent, not a correctness question. What must hold regardless is
            # that dropping (or nearly dropping) it never needs a *larger* ball than Rstar.
            @test cert.support_size <= 4
            @test SB.radius(result) <= Rstar + 1e-7
        elseif delta == 0
            @test SB.radius(result) ≈ Rstar atol = 1e-8
        else
            @test SB.radius(result) >= Rstar
        end
    end
end
