# 1D Hilbert ordering.
#
# In one dimension the Hilbert curve degenerates to the natural interval order: at every
# refinement the left child precedes the right one, so the level-`L` traversal is simply
# `0, 1, ..., 2^L - 1`. A single orientation state therefore suffices, the position of a sector
# is the sector itself, and the state never changes.

const _NSTATES1D = 1

_initialstate(::Val{1}) = 1

_position(::Val{1}, ::Int, sector::Int) = sector

_nextstate(::Val{1}, state::Int, ::Int) = state

_sectorat(::Val{1}, ::Int, position::Int) = position

_nstates(::Val{1}) = _NSTATES1D
