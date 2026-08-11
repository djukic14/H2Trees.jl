# Phase 4: complete 1D randomized oracle.
#
# A 1D ball is an interval [c-r, c+r]. The SEB is [L, U] with L = min(c_i - r_i),
# U = max(c_i + r_i). This is an exact independent oracle for many random instances.

@testset "1D randomized oracle" begin
    rng = MersenneTwister(20240803)
    ntrials = 5000
    maxfail = 0
    for _ in 1:ntrials
        n = rand(rng, 1:8)
        centers = randn(rng, n) .* 10
        radii = abs.(randn(rng, n)) .* 3
        balls = [SB.Ball(SVector(centers[i]), radii[i]) for i in 1:n]
        r = SB.smallest_enclosing_ball(balls)
        cstar, Rstar = oracle_1d(centers, radii)
        scale = max(Rstar, 1.0)
        @test abs(SB.radius(r) - Rstar) <= 1e-7 * scale
        @test abs(SB.center(r)[1] - cstar) <= 1e-7 * scale
        # independent containment
        for i in 1:n
            @test indep_residual(r, balls[i]) <= 1e-7 * scale
        end
    end
    @test maxfail == 0
end

@testset "1D randomized oracle with zero radii mixed in" begin
    rng = MersenneTwister(7)
    for _ in 1:1000
        n = rand(rng, 1:6)
        centers = randn(rng, n) .* 5
        radii = [rand(rng, Bool) ? 0.0 : abs(randn(rng)) for _ in 1:n]
        balls = [SB.Ball(SVector(centers[i]), radii[i]) for i in 1:n]
        r = SB.smallest_enclosing_ball(balls)
        cstar, Rstar = oracle_1d(centers, radii)
        scale = max(Rstar, 1.0)
        @test abs(SB.radius(r) - Rstar) <= 1e-7 * scale
        @test abs(SB.center(r)[1] - cstar) <= 1e-7 * scale
    end
end
