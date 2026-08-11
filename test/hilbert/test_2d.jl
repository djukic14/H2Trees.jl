using Test
using H2Trees

const HO = H2Trees.HilbertOrdering

@testset "2D convention is frozen" begin
    # The documented root convention: sectors visited 0, 1, 3, 2 (the usual U-shape under
    # H2Trees' sector encoding, where dimension 1 is the least significant bit). This is the 2D
    # counterpart of the first x-y portion of the 3D root order.
    @test collect(HO.hilbertsectororder(Val(2), HO.initialstate(Val(2)))) == [0, 1, 3, 2]

    # Level 1 written out as coordinates, so the convention is readable rather than implied by a
    # sector encoding: bottom-left, right, up, left.
    @test [HO.hilbertcoordinates(Val(2), i, 1) for i in 0:3] == [(0, 0), (1, 0), (1, 1), (0, 1)]

    # A full level-2 traversal, frozen explicitly rather than derived at runtime (deriving it
    # from the same tables it is meant to pin would make this test vacuous). Laid out on the
    # grid, with y increasing upward, this is:
    #
    #     15 12 11 10
    #     14 13  8  9
    #      1  2  7  6
    #      0  3  4  5
    @test [HO.hilbertcoordinates(Val(2), i, 2) for i in 0:15] == [
        (0, 0),
        (0, 1),
        (1, 1),
        (1, 0),
        (2, 0),
        (3, 0),
        (3, 1),
        (2, 1),
        (2, 2),
        (3, 2),
        (3, 3),
        (2, 3),
        (1, 3),
        (1, 2),
        (0, 2),
        (0, 3),
    ]
end
