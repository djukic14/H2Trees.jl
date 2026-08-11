using Test
using H2Trees

const HO = H2Trees.HilbertOrdering

@testset "1D is the natural interval order" begin
    @test HO.initialstate(Val(1)) == 1
    @test collect(HO.hilbertsectororder(Val(1), HO.initialstate(Val(1)))) == [0, 1]

    # Single orientation: the state never changes and position == sector at every step.
    for sector in 0:1
        @test HO.hilbertposition(Val(1), 1, sector) == sector
        @test HO.hilbertnextstate(Val(1), 1, sector) == 1
    end

    # The level-L traversal is exactly 0, 1, ..., 2^L - 1, i.e. index == coordinate.
    for level in 1:8
        for x in 0:((1 << level) - 1)
            @test HO.hilbertindex(Val(1), (x,), level) == x
            @test HO.hilbertcoordinates(Val(1), x, level) == (x,)
        end
    end
end
