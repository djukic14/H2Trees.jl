# High-precision (BigFloat) cross-checks for 3D constructed-optimum cases.
#
# This is not an algorithmically independent oracle (both runs use the same SEBB code), but
# running the identical construction at 256-bit precision is useful for catching wrong root
# selection, insufficient discriminant handling, or overly tight/loose affine-independence
# thresholds that only bite at particular scalar-type magnitudes (plan Section 25).

@testset "3D constructed-optimum: Float64 vs BigFloat (well-conditioned tetrahedron)" begin
    setprecision(BigFloat, 256) do
        c64 = SVector(2.0, -3.0, 0.5)
        R64 = 7.0
        dirs64 = tetrahedron_dirs(Float64)
        radii64 = [0.2, 1.1, 2.0, 0.7]
        balls64 = tangent_balls(c64, R64, dirs64, radii64)

        cbig = SVector(BigFloat(2), BigFloat(-3), BigFloat("0.5"))
        Rbig = BigFloat(7)
        dirsbig = tetrahedron_dirs(BigFloat)
        radiibig = BigFloat.(radii64)
        ballsbig = tangent_balls(cbig, Rbig, dirsbig, radiibig)

        r64 = SB.smallest_enclosing_ball(balls64)
        rbig = SB.smallest_enclosing_ball(ballsbig)

        # Each run independently recovers the analytically known (cstar, Rstar).
        @test norm(SB.center(r64) - c64) <= 1e-10
        @test abs(SB.radius(r64) - R64) <= 1e-10
        @test norm(SB.center(rbig) - cbig) <= BigFloat(10)^(-60)
        @test abs(SB.radius(rbig) - Rbig) <= BigFloat(10)^(-60)

        # And the two precisions agree with each other.
        @test norm(SB.center(r64) - Float64.(SB.center(rbig))) <= 1e-9
        @test abs(SB.radius(r64) - Float64(SB.radius(rbig))) <= 1e-9
    end
end

# Squash the (well-conditioned, equal-norm) tetrahedron directions toward the xy-plane by
# scaling their z-component by `epsilon` and renormalizing. Because the four raw directions
# sum to exactly zero component-wise *and* keep an identical norm by symmetry (verified by
# construction: each raw vector is (±1, ±1, ±epsilon) with an even number of sign flips
# canceling in the sum), the squashed, renormalized directions still sum to exactly zero, so
# the equal-radius tangent construction remains the EXACT analytic optimum even as the four
# centers become near-coplanar and the support's Gram matrix becomes near-singular.
function squashed_tetrahedron_dirs(::Type{T}, epsilon) where {T}
    eps_t = T(epsilon)
    raw = SVector{3,T}[
        SVector{3,T}(1, 1, eps_t),
        SVector{3,T}(1, -1, -eps_t),
        SVector{3,T}(-1, 1, -eps_t),
        SVector{3,T}(-1, -1, eps_t),
    ]
    return [d / norm(d) for d in raw]
end

@testset "3D constructed-optimum: Float64 vs BigFloat (near-coplanar support)" begin
    setprecision(BigFloat, 256) do
        cstar64 = SVector(1.0, -1.0, 2.0)
        Rstar64 = 4.0
        radii64 = [0.3, 0.3, 0.3, 0.3]
        epsilon = 1.0e-3

        dirs64 = squashed_tetrahedron_dirs(Float64, epsilon)
        balls64 = tangent_balls(cstar64, Rstar64, dirs64, radii64)

        cstarbig = SVector(BigFloat(1), BigFloat(-1), BigFloat(2))
        Rstarbig = BigFloat(4)
        radiibig = BigFloat.(radii64)
        dirsbig = squashed_tetrahedron_dirs(BigFloat, epsilon)
        ballsbig = tangent_balls(cstarbig, Rstarbig, dirsbig, radiibig)

        r64 = SB.smallest_enclosing_ball(balls64)
        rbig = SB.smallest_enclosing_ball(ballsbig)

        assert_encloses(r64, balls64; rtol=1e-6)
        assert_encloses(rbig, ballsbig; rtol=1e-6)

        # Near-singular support: use a looser tolerance, per plan Section 25.
        @test norm(SB.center(r64) - cstar64) <= 1e-7
        @test abs(SB.radius(r64) - Rstar64) <= 1e-7
        @test norm(SB.center(r64) - Float64.(SB.center(rbig))) <= 1e-7
        @test abs(SB.radius(r64) - Float64(SB.radius(rbig))) <= 1e-7
    end
end
