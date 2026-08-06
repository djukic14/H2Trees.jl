# Called from docs/make.jl before `makedocs`. Regenerates docs/src/assets/plots only if it's
# missing or the content hash of its sources -- PLOT_SOURCE_PATHS (src/, ext/, docs/plots/,
# the driver files, docs/Project.toml) and PLOT_SOURCE_DEPS (resolved versions of the
# external packages the plots depend on), see docs/plotutils.jl -- no longer matches the hash
# recorded the last time it was generated. Nothing under docs/src/assets/plots is committed
# to git, so on a clean checkout this always runs once; on a repeat local build, or in CI when
# the cache step in Documentation.yml hits, it's a no-op.

include(joinpath(@__DIR__, "plotutils.jl"))

const PLOTS_DIR = joinpath(@__DIR__, "src", "assets", "plots")
# Deliberately outside docs/src/assets: that directory is copied verbatim into the deployed
# site, and this hash file is only an internal cache-invalidation marker.
const HASH_FILE = joinpath(@__DIR__, ".plots_hash")

currenthash = sourcehash(joinpath(@__DIR__, ".."))
cachedhash = isfile(HASH_FILE) ? read(HASH_FILE, String) : nothing

if cachedhash == currenthash
    @info "docs/src/assets/plots is up to date, skipping regeneration"
else
    # Evaluated in its own module, not Main: docs/make.jl already has its own `using`
    # statements (e.g. Graphs, for an unrelated example) by the time this runs, and Graphs
    # and CompScienceMeshes both export `vertices` -- running genplots.jl straight into Main
    # would make that ambiguous. genplots.jl only needs what it `using`s itself. Defined via
    # `module ... end` (rather than the bare `Module()` constructor) so it gets its own
    # `include`, which genplots.jl relies on to pull in docs/plotutils.jl and docs/plots/*.jl.
    sandbox = Core.eval(Main, :(module H2TreesDocPlots end))
    Base.include(sandbox, joinpath(@__DIR__, "genplots.jl"))
end
