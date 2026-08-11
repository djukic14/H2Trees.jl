#!/usr/bin/env julia
#
# Pre-renders the interactive PlotlyJS visualizations embedded in the manual pages into
# docs/src/assets/plots/*.html, so `makedocs` no longer has to build meshes/trees and run
# PlotlyJS on every documentation build.
#
# Each figure's code lives in exactly one place, docs/plots/<name>.jl. This script
# `include`s it to produce the plot, and the manual pages `@eval` the same file (via
# `displayedcode` in docs/plotutils.jl) to show its source. There is no second copy of the
# plotting code anywhere, so the two can't drift apart.
#
# docs/src/assets/plots is NOT committed to git: PlotlyJS embeds a fresh random UUID into
# every render, so a committed copy could never be diff-checked for staleness. Normally you
# don't need to run this directly; docs/make.jl calls docs/ensureplots.jl first, which
# only regenerates when the sources covered by `sourcehash` (PLOT_SOURCE_PATHS and
# PLOT_SOURCE_DEPS in docs/plotutils.jl: src/, ext/, docs/plots/, the driver files, and
# resolved dependency versions) actually changed. Run this script directly if you want to
# force a regeneration and look at the result without a full `makedocs` run. See
# docs/src/contributing.md for details.

using CompScienceMeshes
using H2Trees
using PlotlyJS
using ParallelKMeans
using Metis
using BEAST
using Logging

isdefined(@__MODULE__, :sourcehash) || include(joinpath(@__DIR__, "plotutils.jl"))

const PLOTS_SRC = joinpath(@__DIR__, "plots")
const PLOTS_DIR = joinpath(@__DIR__, "src", "assets", "plots")
mkpath(PLOTS_DIR)

# Two of the KMeansTree examples deliberately use an unreliable radius-update policy
# (see the `!!! warning` in h2parallelkmeanstrees.md) which logs a warning on construction;
# that's expected here and not worth surfacing while regenerating plots.
const SUPPRESS_WARNINGS = Set(["kmeans_traceball_2", "kmeans_traceball_3"])

for name in (
    "sphere_tracecube",
    "sphere_traceball",
    "metis_tree_clusters",
    "metis_forest_clusters",
    "kmeans_traceball_1",
    "kmeans_traceball_2",
    "kmeans_traceball_3",
    "tree_families",
    "blocktree",
    "hilbert_curve_1d",
    "hilbert_curve_2d",
    "hilbert_curve_3d",
    "hilbert_tree_1d",
    "hilbert_tree_2d",
    "hilbert_tree_3d",
    "hilbert_tree_adaptive_1d",
    "hilbert_tree_adaptive_2d",
    "hilbert_tree_adaptive_3d",
    "sebb_tangent",
)
    path = joinpath(PLOTS_SRC, "$name.jl")
    p = if name in SUPPRESS_WARNINGS
        with_logger(NullLogger()) do
            return include(path)
        end
    else
        include(path)
    end
    savefig(p, joinpath(PLOTS_DIR, "$name.html"))
end

write(joinpath(@__DIR__, ".plots_hash"), sourcehash(joinpath(@__DIR__, "..")))

@info "Wrote $(length(readdir(PLOTS_DIR))) pre-rendered plot(s) to $PLOTS_DIR"
