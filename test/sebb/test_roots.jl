# `real_roots` is the only scalar equation solver in SEBB, and its degeneracy test decides
# whether a support produces a candidate ball at all. A wrongly-rejected equation is invisible
# from the outside: the enumeration simply runs out of candidates and the whole solve fails.
# These tests therefore check the solver directly against high-precision roots rather than only
# through `smallest_enclosing_ball`.

const ROOTTOL = SB.Tolerance{Float64}(0.0, sqrt(eps(Float64)))

# Independent reference: the same quadratic solved at 256 bits, so cancellation in `Float64`
# cannot hide behind an equally inaccurate expectation.
function reference_roots(A, B, C)
    return setprecision(BigFloat, 256) do
        a, b, c = big(A), big(B), big(C)
        disc = b^2 - 4 * a * c
        disc < 0 && return Float64[]
        s = sqrt(disc)
        return sort([Float64((-b + s) / (2a)), Float64((-b - s) / (2a))])
    end
end

# Relative agreement of a returned root set with the reference.
function roots_match(roots, reference; rtol=1e-10)
    length(reference) == roots.count || return false
    for k in 1:(roots.count)
        scale = max(abs(reference[k]), one(Float64))
        abs(roots.roots[k] - reference[k]) <= rtol * scale || return false
    end
    return true
end

@testset "roots survive a dominant constant term (regression)" begin
    # These are the exact `(A, B, C)` produced by `_candidate_many` for three bounding balls
    # roughly 1e8 units apart: the ordinary situation for a ball tree whose coordinates are
    # not O(1). `A` is O(1) and `C` is O(R^2), so normalizing by the largest coefficient alone
    # made `A` look negligible and the equation was solved as (and then rejected as) a
    # degenerate linear one. No roots came back, no support produced a candidate, and the whole
    # solve failed with "SEBB failed to find a numerically valid enclosing ball".
    A = -0.999999999617245
    B = 993.5408661575898
    C = 1.2499999999506436e16
    roots = SB.real_roots(A, B, C, ROOTTOL)
    @test roots.count == 2
    @test roots_match(roots, reference_roots(A, B, C))
    @test issorted(roots.roots)
end

@testset "roots are recovered across many magnitudes" begin
    # A quadratic with the prescribed roots. The old scaling failed for `|root| > 1/sqrt(rtol)`,
    # i.e. from about 1e4 upward, which is why the exponents here straddle that value.
    for exponent in -8:2:12
        r1 = 0.5 * 10.0^exponent
        r2 = 3.0 * 10.0^exponent
        A = -1.0                       # sign flipped to mirror `_candidate_many`'s convention
        B = r1 + r2
        C = -r1 * r2
        roots = SB.real_roots(A, B, C, ROOTTOL)
        @test roots.count == 2
        @test roots_match(roots, sort([r1, r2]))
    end
end

@testset "scaling every coefficient does not change the roots" begin
    A, B, C = -1.0, 2.5e6, 3.0e11
    base = SB.real_roots(A, B, C, ROOTTOL)
    @test base.count == 2
    for factor in (1.0e-30, 1.0e-8, 7.0, 1.0e12, 1.0e30)
        scaled = SB.real_roots(factor * A, factor * B, factor * C, ROOTTOL)
        @test scaled.count == base.count
        for k in 1:(base.count)
            @test isapprox(scaled.roots[k], base.roots[k]; rtol=1e-12)
        end
    end
end

@testset "degenerate coefficient triples" begin
    # All coefficients zero: every `R` solves it, so no finite root set is meaningful.
    @test SB.real_roots(0.0, 0.0, 0.0, ROOTTOL).count == 0
    # Constant only: unsatisfiable.
    @test SB.real_roots(0.0, 0.0, 5.0, ROOTTOL).count == 0

    # Genuinely linear.
    linear = SB.real_roots(0.0, 4.0, -12.0, ROOTTOL)
    @test linear.count == 1
    @test isapprox(linear.roots[1], 3.0; rtol=1e-12)

    # Linear with a large root, the case the rescaling has to keep exact.
    biglinear = SB.real_roots(0.0, 1.0, -2.5e14, ROOTTOL)
    @test biglinear.count == 1
    @test isapprox(biglinear.roots[1], 2.5e14; rtol=1e-12)

    # Double root: deduplicated to a single entry.
    double = SB.real_roots(1.0, -2.0e6, 1.0e12, ROOTTOL)
    @test double.count == 1
    @test isapprox(double.roots[1], 1.0e6; rtol=1e-7)

    # No real root.
    @test SB.real_roots(1.0, 0.0, 4.0, ROOTTOL).count == 0

    # A linear coefficient small enough that the root is enormous, but still exactly
    # representable: the solver must produce it rather than call the equation degenerate.
    subnormal = SB.real_roots(0.0, 1.0e-300, 1.0, ROOTTOL)
    @test subnormal.count == 1
    @test isapprox(subnormal.roots[1], -1.0e300; rtol=1e-12)
end

@testset "a vanishing quadratic coefficient keeps both roots accurate" begin
    # `A` shrinking towards zero is the one situation a threshold-based "this is really linear"
    # test was meant to protect against. The Citardauq pairing already handles it: one root
    # tends to infinity while the other stays at `-C/B`, and both are computed without
    # cancellation. So the solver must keep returning both, accurately, with no threshold.
    for exponent in (-6, -10, -14, -18)
        A = 10.0^exponent
        B, C = 1.0, -2.0
        roots = SB.real_roots(A, B, C, ROOTTOL)
        reference = reference_roots(A, B, C)
        @test roots.count == 2
        @test roots_match(roots, reference; rtol=1e-9)
        # The small root is the physically meaningful one and must not be swamped by the
        # large. It sits at `-C/B` up to a correction of order `A`, so that is the tolerance.
        @test isapprox(roots.roots[2], 2.0; rtol=10 * A)
    end

    # Exactly zero is the only case that takes the linear branch.
    @test SB.real_roots(0.0, 1.0, -2.0, ROOTTOL).count == 1
    @test isapprox(SB.real_roots(0.0, 1.0, -2.0, ROOTTOL).roots[1], 2.0; rtol=1e-12)

    # A subnormal quadratic coefficient pushes one root past `floatmax`. The overflowing root
    # must be dropped rather than handed on as an `Inf` radius, and the companion root (the
    # one a support could actually use) must survive intact. `2e-310` is chosen so that the
    # normalized coefficient stays *nonzero*: this has to exercise the quadratic branch's
    # overflow handling, not the linear branch.
    @test !iszero(2.0e-310 / 2)
    overflowing = SB.real_roots(2.0e-310, 1.0, -2.0, ROOTTOL)
    @test overflowing.count == 1
    @test all(isfinite, overflowing.roots)
    @test isapprox(overflowing.roots[1], 2.0; rtol=1e-12)

    # One notch smaller and the normalization underflows to exactly zero, which is the linear
    # branch instead. Same answer by a different route, so neither can silently take over.
    @test iszero(5.0e-324 / 2)
    underflowing = SB.real_roots(5.0e-324, 1.0, -2.0, ROOTTOL)
    @test underflowing.count == 1
    @test isapprox(underflowing.roots[1], 2.0; rtol=1e-12)
end

@testset "Float32 roots" begin
    tol32 = SB.Tolerance{Float32}(0.0f0, sqrt(eps(Float32)))
    roots = SB.real_roots(-1.0f0, 4.0f0, -3.0f0, tol32)
    @test roots.count == 2
    @test roots.roots isa SVector{2,Float32}
    @test isapprox(roots.roots[1], 1.0f0; rtol=1.0f-5)
    @test isapprox(roots.roots[2], 3.0f0; rtol=1.0f-5)
end

@testset "a near-zero linear coefficient is not a degenerate equation (regression)" begin
    # The mesh-derived counterpart of the large-coordinate case above: a symmetric support puts
    # `B` at rounding level while `A` and `C` are perfectly ordinary. Any attempt to
    # nondimensionalize `R` by an estimate built from `C / B` blows up here, pushes the roots
    # down to ~1e-15 in the rescaled variable, and the coincident-root test then collapses them
    # to zero: which silently cost the exact answer for 960 real `MetisTree` node balls.
    A = -0.8888888888888891
    B = 3.304678815179371e-17
    C = 0.0012499999999999953
    roots = SB.real_roots(A, B, C, ROOTTOL)
    @test roots.count == 2
    @test roots_match(roots, reference_roots(A, B, C))
    @test isapprox(roots.roots[2], 0.0375; rtol=1e-9)
end
