"""
    HilbertOrdering

Hilbert space-filling-curve helpers used by H2Trees.

The module is tree-independent. It operates only on dimension, orientation state, sector, and
integer grid coordinates.

Supported dimensions are exactly 1, 2, and 3; anything else throws an `ArgumentError` rather
than silently producing an invalid order.

Conventions:

  - sectors are zero-based, `0:2^N-1`, matching `TwoNTree`'s `sectorcentersize` (dimension 1 is
    the least significant bit);
  - Hilbert positions are zero-based, `0:2^N-1`;
  - orientation states are one-based (they index tables) and should be treated as opaque
    outside this module.

Not exported from `H2Trees`; reach it as `H2Trees.HilbertOrdering`.
"""
module HilbertOrdering

using StaticArrays

# Derive traversal order from positions so the inverse table cannot drift.
function _invertpositions(positions::SMatrix{S,NS,Int}) where {S,NS}
    inverted = MMatrix{S,NS,Int}(undef)
    for state in 1:NS
        for sector in 0:(S - 1)
            inverted[positions[sector + 1, state] + 1, state] = sector
        end
    end
    return SMatrix(inverted)
end

include("hilbert1d.jl")
include("hilbert2d.jl")
include("hilbert3d.jl")

# Supported dimensions. Adding one means adding a `hilbert<N>d.jl` with the same four internal
# methods plus a `_checkdimension` method here.
_checkdimension(::Union{Val{1},Val{2},Val{3}}) = nothing

function _checkdimension(::Val{N}) where {N}
    return throw(
        ArgumentError(
            "Hilbert ordering is supported only for dimensions 1, 2, and 3, got $N"
        ),
    )
end

@inline _nsectors(::Val{N}) where {N} = 1 << N

function _checksector(v::Val{N}, sector::Int) where {N}
    0 <= sector < _nsectors(v) || throw(
        ArgumentError(
            "sector must be in 0:$(_nsectors(v) - 1) for dimension $N, got $sector"
        ),
    )
    return nothing
end

function _checkstate(v::Val{N}, state::Int) where {N}
    1 <= state <= _nstates(v) || throw(
        ArgumentError("state must be in 1:$(_nstates(v)) for dimension $N, got $state")
    )
    return nothing
end

function _checklevel(level::Int)
    level >= 0 || throw(ArgumentError("level must be non-negative, got $level"))
    return nothing
end

"""
    initialstate(::Val{N})

Orientation state the Hilbert curve is entered in for dimension `N`.

This is the state that must be handed to the root of a traversal; child states then follow via
[`hilbertnextstate`](@ref).
"""
function initialstate(v::Val{N}) where {N}
    _checkdimension(v)
    return _initialstate(v)
end

"""
    hilbertposition(::Val{N}, state, sector)

Zero-based position at which `sector` is visited by a cell in orientation `state`.

Ordering key for sectors in one cell. For the inverse mapping, use
[`hilbertsectororder`](@ref).
"""
function hilbertposition(v::Val{N}, state::Int, sector::Int) where {N}
    _checkdimension(v)
    _checkstate(v, state)
    _checksector(v, sector)
    return _position(v, state, sector)
end

"""
    hilbertnextstate(::Val{N}, state, sector)

Orientation state of the child in `sector` of a cell in orientation `state`.

Propagating this from the root is what makes the ordering a genuine recursive Hilbert curve
rather than a per-level shuffle.
"""
function hilbertnextstate(v::Val{N}, state::Int, sector::Int) where {N}
    _checkdimension(v)
    _checkstate(v, state)
    _checksector(v, sector)
    return _nextstate(v, state, sector)
end

"""
    hilbertsectororder(::Val{N}, state)

Sectors of a cell in orientation `state`, in Hilbert traversal order.

The inverse permutation of [`hilbertposition`](@ref): entry `k` is the sector visited at
position `k - 1`.
"""
function hilbertsectororder(v::Val{N}, state::Int) where {N}
    _checkdimension(v)
    _checkstate(v, state)
    return ntuple(k -> _sectorat(v, state, k - 1), _nsectors(v))
end

"""
    hilbertindex(::Val{N}, coords, level)

Position of the grid cell at integer `coords` along the level-`level` Hilbert curve.

`coords` holds `N` coordinates, each in `0:2^level-1`. The result is in `0:2^(N*level)-1`.

This is a validation/analysis primitive: tree construction does not need it (it propagates
orientation state directly down the tree instead), and no Hilbert index is ever stored on a
node.
"""
function hilbertindex(v::Val{N}, coords, level::Int) where {N}
    _checkdimension(v)
    _checklevel(level)
    length(coords) == N || throw(
        ArgumentError("expected $N coordinates for dimension $N, got $(length(coords))")
    )
    bound = 1 << level
    for c in coords
        0 <= c < bound || throw(
            ArgumentError("coordinates must be in 0:$(bound - 1) at level $level, got $c"),
        )
    end

    state = _initialstate(v)
    index = 0
    # Refine from the coarsest bit (level 1) down, mirroring a root-to-leaf tree descent: at each
    # step the bit each coordinate contributes forms the sector, exactly as `sectorcentersize`
    # builds one (dimension 1 = least significant bit).
    for step in 1:level
        shift = level - step
        sector = 0
        for d in 1:N
            sector |= ((coords[d] >> shift) & 1) << (d - 1)
        end
        index = (index << N) | _position(v, state, sector)
        state = _nextstate(v, state, sector)
    end
    return index
end

"""
    hilbertcoordinates(::Val{N}, index, level)

Integer grid coordinates of the cell at position `index` along the level-`level` Hilbert curve.

Inverse of [`hilbertindex`](@ref); returns an `NTuple{N,Int}`.
"""
function hilbertcoordinates(v::Val{N}, index::Int, level::Int) where {N}
    _checkdimension(v)
    _checklevel(level)
    ncells = 1 << (N * level)
    0 <= index < ncells ||
        throw(ArgumentError("index must be in 0:$(ncells - 1) at level $level, got $index"))

    state = _initialstate(v)
    coords = zeros(MVector{N,Int})
    mask = _nsectors(v) - 1
    for step in 1:level
        shift = N * (level - step)
        position = (index >> shift) & mask
        sector = _sectorat(v, state, position)
        for d in 1:N
            coords[d] = (coords[d] << 1) | ((sector >> (d - 1)) & 1)
        end
        state = _nextstate(v, state, sector)
    end
    return NTuple{N,Int}(coords)
end

end # module HilbertOrdering
