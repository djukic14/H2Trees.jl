# Phases 2-3: analytically known solutions (one/two/nested/symmetric/equal-radius).

@testset "one ball returns itself" begin
    for c in (SVector(1.0), SVector(1.0, 2.0), SVector(1.0, 2.0, 3.0))
        b = SB.Ball(c, 2.5)
        r = SB.smallest_enclosing_ball([b])
        @test SB.center(r) == c
        @test SB.radius(r) == 2.5
    end
end

@testset "nested balls" begin
    @testset "same center, different radii" begin
        balls = [SB.Ball(SVector(0.0, 0.0), 5.0), SB.Ball(SVector(0.0, 0.0), 1.0)]
        r = SB.smallest_enclosing_ball(balls)
        @test balls_approx(r, balls[1])
    end
    @testset "different centers, one contains all" begin
        big = SB.Ball(SVector(0.0, 0.0, 0.0), 10.0)
        balls = [
            big, SB.Ball(SVector(1.0, 1.0, 1.0), 1.0), SB.Ball(SVector(-2.0, 0.0, 3.0), 0.5)
        ]
        r = SB.smallest_enclosing_ball(balls)
        @test balls_approx(r, big)
        assert_encloses(r, balls)
    end
    @testset "exact internal tangency" begin
        # inner touches outer boundary from inside
        balls = [SB.Ball(SVector(0.0, 0.0), 5.0), SB.Ball(SVector(4.0, 0.0), 1.0)]
        r = SB.smallest_enclosing_ball(balls)
        @test balls_approx(r, balls[1])
    end
    @testset "duplicate largest balls" begin
        big = SB.Ball(SVector(0.0, 0.0), 3.0)
        balls = [big, big, SB.Ball(SVector(0.5, 0.0), 0.5)]
        r = SB.smallest_enclosing_ball(balls)
        @test balls_approx(r, big)
    end
end

@testset "two non-nested balls vs analytic formula" begin
    function two_ball_reference(c1, r1, c2, r2)
        d = norm(c2 .- c1)
        R = (d + r1 + r2) / 2
        t = (R - r1) / d
        return c1 .+ t .* (c2 .- c1), R
    end
    cases = [
        (SVector(0.0), 1.0, SVector(4.0), 1.0),
        (SVector(0.0, 0.0), 1.0, SVector(3.0, 4.0), 2.0),
        (SVector(0.0, 0.0, 0.0), 0.5, SVector(1.0, 2.0, 2.0), 1.5),
        (SVector(-10.0, 5.0), 0.1, SVector(10.0, -5.0), 3.0),
        (SVector(0.0, 0.0), 1.0, SVector(0.001, 0.0), 1.0),  # small separation, non-nested
    ]
    for (c1, r1, c2, r2) in cases
        balls = [SB.Ball(c1, r1), SB.Ball(c2, r2)]
        r = SB.smallest_enclosing_ball(balls)
        cref, Rref = two_ball_reference(c1, r1, c2, r2)
        # For the last (near-nested) case the exact SEB may still be the two-ball ball.
        @test balls_approx(r, SB.Ball(cref, Rref); rtol=1e-6)
        assert_encloses(r, balls)
    end
end

@testset "equal radii == point SEB + r" begin
    # All radii equal r; center is the point-SEB center of the centers, radius = pointR + r.
    pts2 = [SVector(0.0, 0.0), SVector(2.0, 0.0), SVector(1.0, 2.0), SVector(0.5, 0.5)]
    for r0 in (0.0, 0.7, 3.0)
        balls = [SB.Ball(p, r0) for p in pts2]
        pcenter, pradius = BoundingSphere.boundingsphere(pts2)
        res = SB.smallest_enclosing_ball(balls)
        @test balls_approx(res, SB.Ball(SVector(pcenter...), pradius + r0); rtol=1e-6)
        assert_encloses(res, balls)
    end
end

@testset "symmetric: equilateral triangle disks (2D)" begin
    pts = [SVector(cos(2π * k / 3), sin(2π * k / 3)) for k in 0:2]
    for r0 in (0.0, 0.5, 2.0)
        balls = [SB.Ball(p, r0) for p in pts]
        r = SB.smallest_enclosing_ball(balls)
        @test balls_approx(r, SB.Ball(SVector(0.0, 0.0), 1.0 + r0); rtol=1e-6)
        assert_encloses(r, balls)
    end
end

@testset "symmetric: regular tetrahedron (3D)" begin
    tet = [
        SVector(1.0, 1.0, 1.0),
        SVector(1.0, -1.0, -1.0),
        SVector(-1.0, 1.0, -1.0),
        SVector(-1.0, -1.0, 1.0),
    ]
    for r0 in (0.0, 0.3, 1.0)
        balls = [SB.Ball(p, r0) for p in tet]
        r = SB.smallest_enclosing_ball(balls)
        @test balls_approx(r, SB.Ball(SVector(0.0, 0.0, 0.0), sqrt(3.0) + r0); rtol=1e-6)
        assert_encloses(r, balls)
    end
end

@testset "symmetric: square (2D) and cube (3D)" begin
    square = [
        SVector(1.0, 1.0), SVector(1.0, -1.0), SVector(-1.0, 1.0), SVector(-1.0, -1.0)
    ]
    balls = [SB.Ball(p, 0.0) for p in square]
    r = SB.smallest_enclosing_ball(balls)
    @test balls_approx(r, SB.Ball(SVector(0.0, 0.0), sqrt(2.0)); rtol=1e-7)

    cube = [
        SVector(Float64(a), Float64(b), Float64(c)) for a in (-1, 1) for b in (-1, 1) for
        c in (-1, 1)
    ]
    ballsc = [SB.Ball(p, 0.0) for p in cube]
    rc = SB.smallest_enclosing_ball(ballsc)
    @test balls_approx(rc, SB.Ball(SVector(0.0, 0.0, 0.0), sqrt(3.0)); rtol=1e-7)
end

@testset "antipodal pair and extra tangent inputs" begin
    # More tangent inputs than the minimal support size: a diameter defined by two balls with
    # additional balls tangent on the boundary.
    balls = [
        SB.Ball(SVector(-2.0, 0.0), 1.0),
        SB.Ball(SVector(2.0, 0.0), 1.0),
        SB.Ball(SVector(0.0, 3.0), 0.0),   # touches boundary at top of R=3 ball centered origin
        SB.Ball(SVector(0.0, -3.0), 0.0),
    ]
    r = SB.smallest_enclosing_ball(balls)
    @test balls_approx(r, SB.Ball(SVector(0.0, 0.0), 3.0); rtol=1e-7)
    assert_encloses(r, balls)
end
