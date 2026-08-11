# Phase 7: allocation budgets (documented, not zero).
#
# The core enumeration allocates for the canonical-sort copy and small StaticArray work. We
# assert a generous, documented budget rather than requiring zero allocations.

function alloc_of(balls)
    SB.smallest_enclosing_ball(balls)      # warm up
    return @allocated SB.smallest_enclosing_ball(balls)
end

@testset "allocation budget" begin
    b2 = [SB.Ball(SVector(0.0, 0.0, 0.0), 1.0), SB.Ball(SVector(2.0, 0.0, 0.0), 0.5)]
    b4 = [SB.Ball(SVector(Float64(i), Float64(i^2), Float64(-i)), 0.3) for i in 1:4]
    b8 = [SB.Ball(SVector(randn(), randn(), randn()), 0.3) for _ in 1:8]

    a2 = alloc_of(b2)
    a4 = alloc_of(b4)
    a8 = alloc_of(b8)
    @info "SEBB allocations" two = a2 four = a4 eight = a8

    # Budgets: dominated by the canonical Vector{Ball} copy; scale modestly with n.
    @test a2 <= 2_000
    @test a4 <= 4_000
    @test a8 <= 12_000
end
