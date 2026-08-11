using Test
using H2Trees

const HO = H2Trees.HilbertOrdering

# Dimension-generic invariants. These are what establish that each implementation is a genuine
# recursive Hilbert curve rather than merely some self-consistent permutation: the golden tests
# in test_1d/2d/3d.jl pin *which* curve was chosen, these pin that it is a curve at all.
#
# Depths are chosen so the deepest level of each dimension has a few thousand to ~32k cells:
# enough to exercise many orientation transitions while staying a unit test.
const EXHAUSTIVE_LEVELS = Dict(1 => 1:12, 2 => 1:7, 3 => 1:5)
const CONTINUITY_LEVELS = Dict(1 => 1:10, 2 => 1:6, 3 => 1:4)
const INTERVAL_LEVELS = Dict(1 => 1:8, 2 => 1:5, 3 => 1:3)

@testset "state tables are well formed (N=$N)" for N in 1:3
    v = Val(N)
    nsectors = 1 << N
    nstates = N == 1 ? 1 : (N == 2 ? 4 : 12)
    for state in 1:nstates
        positions = [HO.hilbertposition(v, state, sector) for sector in 0:(nsectors - 1)]
        # Positions must be a permutation: two sectors sharing a position would visit one cell
        # twice and skip another.
        @test sort(positions) == collect(0:(nsectors - 1))

        order = collect(HO.hilbertsectororder(v, state))
        @test sort(order) == collect(0:(nsectors - 1))
        # sectororder must be the exact inverse of position, not merely some permutation.
        for sector in 0:(nsectors - 1)
            @test order[HO.hilbertposition(v, state, sector) + 1] == sector
        end

        for sector in 0:(nsectors - 1)
            @test 1 <= HO.hilbertnextstate(v, state, sector) <= nstates
        end
    end
end

@testset "index is a bijection onto 0:2^(N*L)-1 (N=$N)" for N in 1:3
    v = Val(N)
    for level in EXHAUSTIVE_LEVELS[N]
        ncells = 1 << (N * level)
        seen = falses(ncells)
        for index in 0:(ncells - 1)
            coords = HO.hilbertcoordinates(v, index, level)
            recovered = HO.hilbertindex(v, coords, level)
            @test recovered == index
            seen[recovered + 1] = true
        end
        # Detects duplicates and skipped cells, which a per-index round trip alone would not.
        @test all(seen)
    end
end

@testset "consecutive cells share a face (N=$N)" for N in 1:3
    v = Val(N)
    for level in CONTINUITY_LEVELS[N]
        ncells = 1 << (N * level)
        previous = HO.hilbertcoordinates(v, 0, level)
        for index in 1:(ncells - 1)
            current = HO.hilbertcoordinates(v, index, level)
            # L1 distance exactly 1: the spatial-locality property the whole ordering exists for.
            @test sum(abs.(current .- previous)) == 1
            previous = current
        end
    end
end

@testset "children of a cell occupy one contiguous index interval (N=$N)" for N in 1:3
    v = Val(N)
    nsectors = 1 << N
    for level in INTERVAL_LEVELS[N]
        ncells = 1 << (N * level)
        for index in 0:(ncells - 1)
            coords = HO.hilbertcoordinates(v, index, level)
            # The 2^N level-(L+1) cells inside this level-L cell.
            childindices = Int[]
            for sector in 0:(nsectors - 1)
                childcoords = ntuple(d -> (coords[d] << 1) | ((sector >> (d - 1)) & 1), N)
                push!(childindices, HO.hilbertindex(v, childcoords, level + 1))
            end
            # Order within the interval depends on orientation, so compare as a set. This is
            # what distinguishes a genuine recursive curve from one that happens to be
            # continuous at a single level.
            @test Set(childindices) == Set((index * nsectors):((index + 1) * nsectors - 1))
        end
    end
end

@testset "invalid input is rejected" begin
    for N in (0, 4, -1)
        @test_throws ArgumentError HO.initialstate(Val(N))
        @test_throws ArgumentError HO.hilbertposition(Val(N), 1, 0)
        @test_throws ArgumentError HO.hilbertnextstate(Val(N), 1, 0)
        @test_throws ArgumentError HO.hilbertsectororder(Val(N), 1)
        @test_throws ArgumentError HO.hilbertindex(Val(N), ntuple(_ -> 0, max(N, 0)), 1)
        @test_throws ArgumentError HO.hilbertcoordinates(Val(N), 0, 1)
    end

    @test_throws ArgumentError HO.hilbertposition(Val(3), 1, -1)
    @test_throws ArgumentError HO.hilbertposition(Val(3), 1, 8)
    @test_throws ArgumentError HO.hilbertposition(Val(2), 1, 4)
    @test_throws ArgumentError HO.hilbertposition(Val(3), 0, 0)
    @test_throws ArgumentError HO.hilbertposition(Val(3), 13, 0)
    @test_throws ArgumentError HO.hilbertnextstate(Val(2), 5, 0)

    @test_throws ArgumentError HO.hilbertindex(Val(3), (0, 0, 0), -1)
    # coordinate outside 0:2^L-1
    @test_throws ArgumentError HO.hilbertindex(Val(3), (0, 0, 4), 2)
    @test_throws ArgumentError HO.hilbertindex(Val(3), (-1, 0, 0), 2)
    # wrong number of coordinates
    @test_throws ArgumentError HO.hilbertindex(Val(3), (0, 0), 2)

    @test_throws ArgumentError HO.hilbertcoordinates(Val(3), -1, 2)
    @test_throws ArgumentError HO.hilbertcoordinates(Val(3), 1 << 6, 2)
end

@testset "primitives are allocation free" begin
    # These run once per child per node during construction, so they must stay inlineable and
    # non-allocating; a table accidentally becoming a Vector-of-Vector would show up here.
    for N in 1:3
        v = Val(N)
        state = HO.initialstate(v)
        HO.hilbertposition(v, state, 0)
        HO.hilbertnextstate(v, state, 0)
        @test @allocated(HO.hilbertposition(v, state, 0)) == 0
        @test @allocated(HO.hilbertnextstate(v, state, 0)) == 0
    end
end
