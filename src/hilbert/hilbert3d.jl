# 3D Hilbert ordering.
#
# These are the twelve-state tables H2Trees has always used, moved here verbatim from
# `src/trees/clustertrees.jl` (where they lived as `hilbert_positions`/`hilbert_states`). They
# are NOT regenerated from a generic construction: 3D Hilbert curves are not unique, so a
# mathematically valid generic implementation would produce a different, equally valid ordering
# and would silently change every existing H2Trees tree layout. `test/hilbert/test_3d.jl` holds
# a frozen copy of the original tables and asserts this module reproduces them entry for entry.
#
# The only change from the original is that next-state values are stored 1-based here. The old
# tables stored them 0-based and every call site wrote `hilbert_states[state][sector + 1] + 1`;
# folding that `+ 1` into the table keeps the offset in one place instead of at each use.
#
# As in 2D, `_POSITIONS3D[sector + 1, state]` is the position at which a sector is visited, and
# the traversal ORDER is the inverse permutation. Here the two genuinely differ: the root state's
# positions are `0, 1, 3, 2, 7, 6, 4, 5` while its sector order is `0, 1, 3, 2, 6, 7, 5, 4`.

const _NSTATES3D = 12

# Column = orientation state, row = sector + 1.
const _POSITIONS3D = SMatrix{8,12,Int}(
    #= state  1 =#
    0,
    1,
    3,
    2,
    7,
    6,
    4,
    5,
    #= state  2 =#
    0,
    7,
    1,
    6,
    3,
    4,
    2,
    5,
    #= state  3 =#
    0,
    3,
    7,
    4,
    1,
    2,
    6,
    5,
    #= state  4 =#
    2,
    3,
    1,
    0,
    5,
    4,
    6,
    7,
    #= state  5 =#
    4,
    3,
    5,
    2,
    7,
    0,
    6,
    1,
    #= state  6 =#
    6,
    5,
    1,
    2,
    7,
    4,
    0,
    3,
    #= state  7 =#
    4,
    7,
    3,
    0,
    5,
    6,
    2,
    1,
    #= state  8 =#
    6,
    7,
    5,
    4,
    1,
    0,
    2,
    3,
    #= state  9 =#
    2,
    5,
    3,
    4,
    1,
    6,
    0,
    7,
    #= state 10 =#
    2,
    1,
    5,
    6,
    3,
    0,
    4,
    7,
    #= state 11 =#
    4,
    5,
    7,
    6,
    3,
    2,
    0,
    1,
    #= state 12 =#
    6,
    1,
    7,
    0,
    5,
    2,
    4,
    3,
)

const _NEXTSTATES3D = SMatrix{8,12,Int}(
    #= state  1 =#
    2,
    3,
    4,
    3,
    5,
    6,
    4,
    6,
    #= state  2 =#
    3,
    7,
    1,
    8,
    9,
    9,
    1,
    8,
    #= state  3 =#
    1,
    10,
    11,
    10,
    2,
    2,
    12,
    12,
    #= state  4 =#
    7,
    1,
    7,
    12,
    10,
    1,
    10,
    9,
    #= state  5 =#
    12,
    12,
    1,
    8,
    6,
    10,
    1,
    8,
    #= state  6 =#
    5,
    5,
    9,
    9,
    1,
    7,
    11,
    7,
    #= state  7 =#
    6,
    8,
    6,
    4,
    2,
    2,
    12,
    12,
    #= state  8 =#
    7,
    2,
    7,
    11,
    10,
    5,
    10,
    11,
    #= state  9 =#
    11,
    4,
    2,
    2,
    11,
    4,
    6,
    10,
    #= state 10 =#
    5,
    5,
    9,
    9,
    3,
    8,
    3,
    4,
    #= state 11 =#
    8,
    3,
    12,
    3,
    8,
    6,
    9,
    6,
    #= state 12 =#
    11,
    4,
    3,
    7,
    11,
    4,
    5,
    5,
)

const _SECTORORDER3D = _invertpositions(_POSITIONS3D)

_initialstate(::Val{3}) = 1

_position(::Val{3}, state::Int, sector::Int) = _POSITIONS3D[sector + 1, state]

_nextstate(::Val{3}, state::Int, sector::Int) = _NEXTSTATES3D[sector + 1, state]

_sectorat(::Val{3}, state::Int, position::Int) = _SECTORORDER3D[position + 1, state]

_nstates(::Val{3}) = _NSTATES3D
