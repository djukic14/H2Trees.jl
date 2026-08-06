# Shared by docs/genplots.jl (which executes docs/plots/*.jl to produce the pre-rendered
# assets in docs/src/assets/plots/) and the manual pages (which use `displayedcode` in an
# `@eval` block to show the exact same source verbatim). Single source of truth per figure:
# there is no separate copy of the plotting code living in the .md files.

using Markdown
using Pkg

# Text between these two markers still runs when the file is `include`d, but is not shown on
# the manual page -- it's rendering boilerplate (trace/layout construction), or setup that's
# already shown verbatim earlier on the same page, and would just distract from whatever the
# page is actually teaching (Metis partitioning, k-means, tree families, ...). If only the
# start marker is present, everything from it to the end of the file is hidden.
const DOCS_HIDE_START = "# --- hide-from-docs ---"
const DOCS_HIDE_END = "# --- end-hide-from-docs ---"

function displayedcode(path::AbstractString)
    text = read(path, String)
    startidx = findfirst(DOCS_HIDE_START, text)
    if isnothing(startidx)
        shown = text
    else
        before = text[1:(first(startidx) - 1)]
        endidx = findfirst(DOCS_HIDE_END, text)
        after = isnothing(endidx) ? "" : text[(last(endidx) + 1):end]
        shown = isempty(strip(after)) ? before : strip(before) * "\n\n" * strip(after)
    end
    return Markdown.parse("```julia\n" * strip(shown) * "\n```")
end

# The pre-rendered plots in docs/src/assets/plots are not committed to git (PlotlyJS embeds a
# fresh random UUID into every render, so committed HTML can never be diffed for staleness --
# see the "Documentation plots" section of contributing.md). Instead they're reconstructed on
# demand, guarded by a hash over everything that can change what they look like: the plot
# scripts themselves, all of H2Trees (core + extensions, since the scripts build/iterate real
# trees), and the driver/docs-env files that control how they're generated. Each entry is
# either a directory (hashed recursively) or a single file.
const PLOT_SOURCE_PATHS = [
    "src",
    "ext",
    "docs/plots",
    "docs/plotutils.jl",
    "docs/genplots.jl",
    "docs/ensureplots.jl",
    "docs/Project.toml",
]

# External packages the plot scripts depend on but whose source isn't in this repo, so
# wouldn't otherwise be covered by PLOT_SOURCE_PATHS -- a version bump (e.g. a PlotlyJS
# release changing its HTML/JS output) should also invalidate the cache.
const PLOT_SOURCE_DEPS = [
    "PlotlyJS",
    "PlotlyBase",
    "CompScienceMeshes",
    "BEAST",
    "Metis",
    "ParallelKMeans",
    "H2Trees",
]

function depsignature()
    io = IOBuffer()
    for info in sort(collect(values(Pkg.dependencies())); by=d -> d.name)
        info.name in PLOT_SOURCE_DEPS || continue
        print(io, info.name, "@", something(info.version, "dev"), ";")
    end
    return String(take!(io))
end

function sourcehash(pkgroot::AbstractString)
    h = zero(UInt64)
    for entry in PLOT_SOURCE_PATHS
        fullpath = joinpath(pkgroot, entry)
        if isfile(fullpath)
            h = hash(read(fullpath), h)
            h = hash(relpath(fullpath, pkgroot), h)
        elseif isdir(fullpath)
            for (dir, _, files) in walkdir(fullpath)
                for file in sort(files)
                    filepath = joinpath(dir, file)
                    h = hash(read(filepath), h)
                    h = hash(relpath(filepath, pkgroot), h)
                end
            end
        end
    end
    h = hash(depsignature(), h)
    return string(h; base=16)
end
