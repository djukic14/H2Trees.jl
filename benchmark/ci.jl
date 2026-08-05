# Compares `benchmark/benchmarks.jl`'s SUITE between two git revisions of this package: the PR
# head SHA (`target`) and the PR base SHA (`baseline`), read from `ARGS` so
# `.github/workflows/Benchmarks.yml` can pass the two commits a `pull_request` event actually
# diffs. Falls back to `HEAD` vs `origin/main` so this also works for a manual `workflow_dispatch`
# run with no PR context.
#
# Both revisions run consecutively in this same process/machine (PkgBenchmark checks each one out
# of this same repo in turn, benchmarks it, then restores the original checkout) -- deliberately
# not comparing against a report generated on a different day, since runner hardware, OS image,
# and machine load can differ between separate CI runs.
#
# `baseline` might not have `benchmark/benchmarks.jl` yet (the PR that first adds this suite has
# no such baseline to compare against) -- detected below via `hasbenchmarksuite` and handled by
# reporting a target-only run instead of failing. If `target` itself lacks the suite (e.g. a PR
# that removes it), that is a real problem and is left to fail naturally: `judge`/`benchmarkpkg`
# resolve `benchmark/benchmarks.jl` from the currently checked-out tree, which is `target`'s.

using Pkg
Pkg.activate(@__DIR__)

using PkgBenchmark

const PKGDIR = dirname(@__DIR__)

target = length(ARGS) >= 1 ? ARGS[1] : "HEAD"
baseline = length(ARGS) >= 2 ? ARGS[2] : "origin/main"

# One thread: more reproducible than a multi-threaded comparison, and makes an algorithmic
# regression easier to tell apart from scheduling/contention noise (see the benchmark plan's
# "Threaded performance" section -- a compact multi-threaded comparison is a deliberate,
# separate follow-up, not part of this initial suite).
const BENCH_ENV = Dict("JULIA_NUM_THREADS" => "1")

"""
    hasbenchmarksuite(pkgdir, rev)

Whether `benchmark/benchmarks.jl` exists at git revision `rev` of the repo at `pkgdir`, without
checking anything out.
"""
function hasbenchmarksuite(pkgdir, rev)
    return success(
        pipeline(
            `git -C $pkgdir show $(rev):benchmark/benchmarks.jl`;
            stdout=devnull,
            stderr=devnull,
        ),
    )
end

reportfile = joinpath(PKGDIR, "benchmark-report.md")

println("Comparing target=", target, " against baseline=", baseline)

if hasbenchmarksuite(PKGDIR, baseline)
    judgement = judge(
        PKGDIR,
        BenchmarkConfig(; id=target, env=BENCH_ENV),
        BenchmarkConfig(; id=baseline, env=BENCH_ENV);
        judgekwargs=Dict(:time_tolerance => 0.10, :memory_tolerance => 0.01),
    )
    # `export_invariants=true` so the report also lists cases that did NOT regress/improve -- an
    # all-green comparison would otherwise export an empty file, which is not a useful job summary.
    export_markdown(reportfile, judgement; export_invariants=true)
else
    println(
        "Baseline ",
        baseline,
        " has no benchmark/benchmarks.jl -- reporting target-only results instead of a",
        " comparison.",
    )
    results = benchmarkpkg(PKGDIR, BenchmarkConfig(; id=target, env=BENCH_ENV))
    open(reportfile, "w") do io
        println(io, "# Benchmark Report (target-only)")
        println(io)
        println(
            io,
            "No comparison was run: the baseline revision (`",
            baseline,
            "`) does not have `benchmark/benchmarks.jl` yet. Once it merges, later PRs will get",
            " a normal target-vs-baseline comparison.",
        )
        println(io)
        return export_markdown(io, results)
    end
end

println("Wrote ", reportfile)
