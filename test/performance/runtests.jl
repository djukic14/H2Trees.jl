module TestPerformanceContracts

# Performance-contract suite: inference and allocation budgets for tree construction, iterators,
# plans, and `checkadmissibility`.
#
# These complement the correctness tests -- they check HOW the result was computed (concrete
# types, bounded allocations), not whether it is correct. Compiler and allocation behavior is more
# version-sensitive than correctness, so this suite is kept separate from `test/runtests.jl`'s
# ordinary tests (see the "Performance contracts" group there) and is runnable on its own:
#
#     julia --project=. -e 'include("test/performance/runtests.jl")'
#
# No wall-clock timing is asserted anywhere in this suite -- only inference and allocation
# counts/bytes.
#
# Exact allocation counts/bytes can shift between Julia minor versions (compiler changes)
# independently of anything in this package -- the strict allocation checks (`allocations.jl`) are
# therefore pinned to the Julia minor version they were measured and tuned against, matching CI's
# primary version. `inference.jl` (`@inferred`) is not pinned: concrete-type inference is a
# language-semantics guarantee, not a compiler-heuristic/allocator detail, so it is far less likely
# to drift between minor versions.
const PERF_PINNED_JULIA_MINOR_VERSION = v"1.12"
const PERF_STRICT_CHECKS_ACTIVE =
    VERSION.major == PERF_PINNED_JULIA_MINOR_VERSION.major &&
    VERSION.minor == PERF_PINNED_JULIA_MINOR_VERSION.minor

using Test

include("fixtures.jl")
include("workloads.jl")
include("budgets.jl")

@testset verbose = true "Performance contracts" begin
    @testset "Inference" begin
        include("inference.jl")
    end
    if PERF_STRICT_CHECKS_ACTIVE
        @testset "Allocations" begin
            include("allocations.jl")
        end
    else
        @info "Skipping allocation-budget checks: pinned to Julia $(PERF_PINNED_JULIA_MINOR_VERSION.major).$(PERF_PINNED_JULIA_MINOR_VERSION.minor), running on $VERSION"
    end
end

end # module TestPerformanceContracts
