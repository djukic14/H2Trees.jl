# 2D Hilbert ordering.
#
# Four orientation states. The root state (2) traverses sectors in the order
#
#     0, 1, 3, 2
#
# which, with H2Trees' sector convention (dimension 1 is the least significant bit, see
# `sectorcentersize`), is the usual U-shaped Hilbert traversal. That convention is deliberately
# frozen: it is the 2D counterpart of the first x-y portion of the existing 3D root order
# (`0, 1, 3, 2, 6, 7, 5, 4`), and `test/hilbert/test_2d.jl` pins it so a future edit cannot
# silently swap in a different (also valid) 2D Hilbert curve.
#
# NOTE the distinction between the two tables' meaning, which is easy to conflate:
#   * `_POSITIONS2D[sector + 1, state]` is the POSITION AT WHICH that sector is visited;
#   * the traversal ORDER (which sector comes first, second, ...) is the INVERSE permutation,
#     `_SECTORORDER2D`, derived below rather than written out by hand.
# For the 2D root state the two happen to coincide (the permutation is an involution); in 3D
# they genuinely differ, so nothing here may rely on them being interchangeable.

const _NSTATES2D = 4

# Column = orientation state, row = sector + 1.
const _POSITIONS2D = SMatrix{4,4,Int}(
    #= state 1 =#
    0,
    3,
    1,
    2,
    #= state 2 =#
    0,
    1,
    3,
    2,
    #= state 3 =#
    2,
    3,
    1,
    0,
    #= state 4 =#
    2,
    1,
    3,
    0,
)

const _NEXTSTATES2D = SMatrix{4,4,Int}(
    #= state 1 =#
    2,
    3,
    1,
    1,
    #= state 2 =#
    1,
    2,
    4,
    2,
    #= state 3 =#
    3,
    1,
    3,
    4,
    #= state 4 =#
    4,
    4,
    2,
    3,
)

const _SECTORORDER2D = _invertpositions(_POSITIONS2D)

# State 2, not 1: state 1 is the table's first row but the curve is entered in the orientation
# whose sector order is `0, 1, 3, 2` (see above).
_initialstate(::Val{2}) = 2

_position(::Val{2}, state::Int, sector::Int) = _POSITIONS2D[sector + 1, state]

_nextstate(::Val{2}, state::Int, sector::Int) = _NEXTSTATES2D[sector + 1, state]

_sectorat(::Val{2}, state::Int, position::Int) = _SECTORORDER2D[position + 1, state]

_nstates(::Val{2}) = _NSTATES2D
