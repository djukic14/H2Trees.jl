
# Contributing

In order to contribute to this package directly create a pull request against the `main` branch. Before doing so please:  

- Follow the style of the surrounding code.
- Supplement the documentation.
- Write tests and check that no errors occur.

---

## Style

For a consistent style the [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) package is used which enforces the style defined in the *.JuliaFormatter.toml* file. To follow this style simply run

```julia
using JuliaFormatter
format(pkgdir(H2Trees; overwrite=true))
```

!!! note
    That all files follow the JuliaFormatter style is tested during the unit tests. Hence, do not forget to execute the two lines above. Otherwise, the tests are likely to not pass.

---

## Documentation

Add documentation for any changes or new features following the style of the existing documentation. For more information you can have a look at the [Documenter.jl](https://documenter.juliadocs.org/stable/) documentation.

---

## Tests

Write tests for your code changes and verify that no errors occur, e.g., by running

```julia
using Pkg
Pkg.test("H2Trees")
```

For more detailed information on which parts are tested the coverage can be evaluated on your local machine, e.g., by

```julia
using Pkg
Pkg.test("H2Trees"; coverage=true, julia_args=`--threads 6`)

# determine coverage
using Coverage
src_folder = pkgdir(H2Trees) * "/src"
coverage   = process_folder(src_folder)
LCOV.writefile("path-to-folder-you-like" * "H2Trees.lcov.info", coverage)

clean_folder(src_folder) # delete .cov files

# extract information about coverage
covered_lines, total_lines = get_summary(coverage)
@info "Current coverage:\n$covered_lines of $total_lines lines ($(round(Int, covered_lines / total_lines * 100)) %)"
```

In Visual Studio Code the [Coverage Gutters](https://marketplace.visualstudio.com/items?itemName=ryanluker.vscode-coverage-gutters) plugin can be used to visualize the tested lines of the code by inserting the path of the *H2Trees.lcov.info* file in the settings.

---

## Documentation plots

The interactive PlotlyJS visualizations in the manual (`tree_families.md`, `blocktree_petrov.md`,
and the `ext/h2plotlyjstrees.md`/`h2metistrees.md`/`h2parallelkmeanstrees.md` extension pages) are
pre-rendered into `docs/src/assets/plots/*.html` instead of being rebuilt by `makedocs` every time.

Each figure's code lives in exactly one place, `docs/plots/<name>.jl` — there is no copy of it in
the `.md` files. `docs/genplots.jl` `include`s each of these to produce the plot, and the manual
pages use a `@eval` block (via `displayedcode` in `docs/plotutils.jl`) to show that same file's
source verbatim, so the displayed code and the code that actually produced the picture cannot
drift apart. A `# --- hide-from-docs ---` / `# --- end-hide-from-docs ---` marker pair in a
`docs/plots/*.jl` file excludes the rendering boilerplate (or setup already shown earlier on the
page) from what's displayed, without excluding it from execution.

`docs/src/assets/plots/` is **not** committed to git: `PlotlyJS.Plot()` embeds a fresh random UUID
into every render, so a committed copy could never be diffed for staleness — any regeneration
looks "changed" even with byte-identical inputs. Instead, `docs/make.jl` calls
`docs/ensureplots.jl` before `makedocs`, which regenerates only if the assets are missing or the
hash from `sourcehash` (in `docs/plotutils.jl`) no longer matches the one recorded the last time
it ran — so a normal `julia --project=docs docs/make.jl` always produces up-to-date pictures with
no extra step on your part. That hash covers everything that can change what a plot looks like:
all of `src/` and `ext/` (the plot scripts build and iterate real `H2Trees` trees, not just call
the plotting extensions), `docs/plots/`, the driver files themselves (`genplots.jl`,
`plotutils.jl`, `ensureplots.jl`, `docs/Project.toml`), and the resolved versions of the external
packages the plots depend on (`PlotlyJS`, `PlotlyBase`, `CompScienceMeshes`, `BEAST`, `Metis`,
`ParallelKMeans`) — so e.g. a `PlotlyJS` release changing its HTML output also invalidates the
cache, even though its source isn't in this repo. In CI, `Documentation.yml` additionally caches
`docs/src/assets/plots` keyed on the same file set (`hashFiles(...)` can't see resolved dependency
versions, only file contents) so unrelated PRs skip the ~1 minute regeneration cost entirely; if
that cache key happens to hit while the *resolved* dependency versions have actually moved on, the
`sourcehash` check inside `ensureplots.jl` still catches it and regenerates anyway.

If you want to inspect a regenerated plot directly, without a full `makedocs` run:

```julia
julia --project=docs docs/genplots.jl
```
