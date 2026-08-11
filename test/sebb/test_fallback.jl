# The approximate fallback.
#
# `smallest_enclosing_ball` must never fail on a valid input: it runs deep inside tree
# construction, where an exception destroys a whole build over a rounding accident. When the
# exact enumeration cannot resolve a configuration it warns and returns a ball that is
# guaranteed to enclose but is not certified minimal. These tests pin both halves of that
# contract: the guarantee, and the fact that the ordinary path does NOT take it.

using Logging

# Balls the exact enumeration genuinely cannot resolve at `Float64`: three points sitting a
# handful of ULPs apart around 1e12, so the whole geometry is narrower than `eps` at that
# magnitude. `norm(c1 - c2)` then carries an absolute error comparable to the answer itself,
# every tangency residual is noise, and no candidate can pass a containment check. The inputs
# have already lost the information the exact solver would need; nothing downstream can recover
# it. Found by fuzzing; kept verbatim so the fallback path stays reachable from the suite.
const UNRESOLVABLE = [
    SB.Ball(SVector(9.999999999999982e11), 0.0),
    SB.Ball(SVector(9.999999999999989e11), 0.0),
    SB.Ball(SVector(1.0000000000000002e12), 0.0),
]

@testset "the fallback is reachable and warns" begin
    # `maxlog=1` keeps a failing tree build from emitting one warning per node. The limit is
    # counted per logger, so a freshly built `TestLogger` always sees the first message and
    # this does not depend on whether an earlier test already triggered the fallback.
    logger = Test.TestLogger()
    cert = with_logger(logger) do
        return SB._smallest_enclosing_ball_with_certificate(UNRESOLVABLE)
    end
    @test any(r -> r.level == Logging.Warn && occursin("SEBB", r.message), logger.logs)

    # A certificate-free result is how the fallback reports itself.
    @test !SB.isexactresult(cert)
    @test cert.support_size == 0
    @test isempty(SB.active_support(cert))
    assert_encloses(cert.ball, UNRESOLVABLE; rtol=0, atol=0)
end

@testset "the ordinary path does not take the fallback" begin
    # Guards against the fallback quietly becoming the normal answer: every one of these is a
    # configuration the exact enumeration is expected to solve, so all must carry a support.
    logger = Test.TestLogger()
    rng = MersenneTwister(4242)
    with_logger(logger) do
        for N in 1:3, _ in 1:400
            n = rand(rng, 1:6)
            scale = 10.0^rand(rng, -4:8)
            balls = [
                SB.Ball(
                    SVector{N,Float64}(scale .* randn(rng, N)), scale * abs(randn(rng))
                ) for _ in 1:n
            ]
            cert = SB._smallest_enclosing_ball_with_certificate(balls)
            @test SB.isexactresult(cert)
            @test 1 <= cert.support_size <= N + 1
        end
    end
    @test isempty(logger.logs)
end

@testset "_fallbackball always encloses and is close to optimal" begin
    rng = MersenneTwister(31415)
    worstratio = 1.0
    for N in 1:3, _ in 1:300
        n = rand(rng, 1:8)
        balls = [SB.Ball(SVector{N,Float64}(randn(rng, N)), abs(randn(rng))) for _ in 1:n]
        approximate = SB._fallbackball(balls)
        # Containment is a guarantee, not a tolerance: the radius IS the measured enclosing
        # radius at that center, so a strict comparison must hold.
        for b in balls
            @test indep_residual(approximate, b) <= 0
        end
        exact = SB.smallest_enclosing_ball(balls)
        @test SB.radius(approximate) >= SB.radius(exact) * (1 - 1e-9)
        worstratio = max(worstratio, SB.radius(approximate) / max(SB.radius(exact), eps()))
    end
    # Documented quality of the Bădoiu-Clarkson refinement at `_FALLBACKSTEPS` steps.
    @test worstratio <= 1.01
end

@testset "_fallbackball degeneracies" begin
    # Single ball: the answer is that ball, and the farthest-ball direction is undefined
    # (the iterate sits exactly on its center), which must not produce a NaN.
    one = [SB.Ball(SVector(1.0, -2.0, 0.5), 3.0)]
    @test balls_approx(SB._fallbackball(one), one[1])

    # Identical balls: same answer, still finite.
    same = fill(SB.Ball(SVector(0.0, 0.0), 2.0), 4)
    @test balls_approx(SB._fallbackball(same), same[1])

    # Concentric, nested: the containing ball is the answer, and the two-ball seed has to clamp
    # its tangent point back inside the segment to find it.
    nested = [SB.Ball(SVector(0.0, 0.0), 5.0), SB.Ball(SVector(0.5, 0.0), 0.25)]
    approximate = SB._fallbackball(nested)
    @test balls_approx(approximate, nested[1]; rtol=1e-6)

    # Coincident centers, different radii: `_diametralcenter`'s zero-distance branch.
    coincident = [SB.Ball(SVector(3.0, 3.0), 1.0), SB.Ball(SVector(3.0, 3.0), 4.0)]
    @test balls_approx(SB._fallbackball(coincident), coincident[2]; rtol=1e-9)

    # Zero-radius points only: reduces to the smallest enclosing ball of a point set.
    points = [SB.Ball(SVector(-1.0, 0.0), 0.0), SB.Ball(SVector(1.0, 0.0), 0.0)]
    @test balls_approx(SB._fallbackball(points), SB.Ball(SVector(0.0, 0.0), 1.0); rtol=1e-6)
end

@testset "the public API does not throw on any valid input" begin
    # The contract the tree builders rely on. Only genuinely malformed input (no balls,
    # mismatched lengths) may raise.
    rng = MersenneTwister(2718)
    logger = Test.TestLogger()
    with_logger(logger) do
        for _ in 1:2000
            N = rand(rng, 1:3)
            n = rand(rng, 1:8)
            # Deliberately nasty: exact lattice coordinates over a huge dynamic range, with
            # many zero radii, which is where the enumeration is most likely to give up.
            scale = 10.0^rand(rng, -10:10)
            balls = [
                SB.Ball(
                    SVector{N,Float64}(ntuple(_ -> scale * rand(rng, -1:1), N)),
                    rand(rng, Bool) ? 0.0 : scale * rand(rng) * 1.0e-3,
                ) for _ in 1:n
            ]
            result = SB.smallest_enclosing_ball(balls)
            @test result isa SB.Ball{N,Float64}
            for b in balls
                @test indep_residual(result, b) <=
                    sqrt(eps(Float64)) * max(SB.radius(result), 1.0)
            end
        end
    end
end
