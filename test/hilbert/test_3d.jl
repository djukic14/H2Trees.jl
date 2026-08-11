using Test
using H2Trees

const HO = H2Trees.HilbertOrdering

# Frozen copy of the tables that lived in `src/trees/clustertrees.jl` before the Hilbert
# extraction, verbatim and still 0-based in the state column (the module stores next states
# 1-based; the conversion is applied where they are compared below).
#
# The point of freezing them here is that 3D Hilbert curves are NOT unique: a generic
# implementation could produce a different, perfectly valid 3D ordering and silently relabel
# every existing H2Trees tree. These tests fail if that ever happens.
const FROZEN_STATES_3D = [
    [1, 2, 3, 2, 4, 5, 3, 5],
    [2, 6, 0, 7, 8, 8, 0, 7],
    [0, 9, 10, 9, 1, 1, 11, 11],
    [6, 0, 6, 11, 9, 0, 9, 8],
    [11, 11, 0, 7, 5, 9, 0, 7],
    [4, 4, 8, 8, 0, 6, 10, 6],
    [5, 7, 5, 3, 1, 1, 11, 11],
    [6, 1, 6, 10, 9, 4, 9, 10],
    [10, 3, 1, 1, 10, 3, 5, 9],
    [4, 4, 8, 8, 2, 7, 2, 3],
    [7, 2, 11, 2, 7, 5, 8, 5],
    [10, 3, 2, 6, 10, 3, 4, 4],
]

const FROZEN_POSITIONS_3D = [
    [0, 1, 3, 2, 7, 6, 4, 5],
    [0, 7, 1, 6, 3, 4, 2, 5],
    [0, 3, 7, 4, 1, 2, 6, 5],
    [2, 3, 1, 0, 5, 4, 6, 7],
    [4, 3, 5, 2, 7, 0, 6, 1],
    [6, 5, 1, 2, 7, 4, 0, 3],
    [4, 7, 3, 0, 5, 6, 2, 1],
    [6, 7, 5, 4, 1, 0, 2, 3],
    [2, 5, 3, 4, 1, 6, 0, 7],
    [2, 1, 5, 6, 3, 0, 4, 7],
    [4, 5, 7, 6, 3, 2, 0, 1],
    [6, 1, 7, 0, 5, 2, 4, 3],
]

@testset "3D regression against the pre-refactor tables" begin
    # Every (state, sector) pair must resolve exactly as the old tables did. The old call sites
    # wrote `hilbert_states[state][sector + 1] + 1`, so the `+ 1` appears on the frozen side.
    for state in 1:12, sector in 0:7
        @test HO.hilbertposition(Val(3), state, sector) ==
            FROZEN_POSITIONS_3D[state][sector + 1]
        @test HO.hilbertnextstate(Val(3), state, sector) ==
            FROZEN_STATES_3D[state][sector + 1] + 1
    end

    @test HO.initialstate(Val(3)) == 1

    # The documented root convention. NOTE this is the traversal ORDER, i.e. the inverse of the
    # positions row `[0, 1, 3, 2, 7, 6, 4, 5]`. The two differ in 3D, and conflating
    # them is the easiest way to silently install a wrong curve.
    @test collect(HO.hilbertsectororder(Val(3), HO.initialstate(Val(3)))) ==
        [0, 1, 3, 2, 6, 7, 5, 4]
end
