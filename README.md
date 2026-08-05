# 𝓗² Trees

<p align="center">
<picture>
  <source media="(prefers-color-scheme)" srcset="docs/src/assets/logo.svg" height="190">
  <img alt="" src="" height="190">
</picture>
</p>

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://djukic14.github.io/H2Trees.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://djukic14.github.io/H2Trees.jl/dev/)
[![Build Status](https://github.com/djukic14/H2Trees.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/djukic14/H2Trees.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/djukic14/H2Trees.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/djukic14/H2Trees.jl)

## Introduction

The [H2Trees](https://github.com/djukic14/H2Trees.jl) package provides Julia tree data
structures, traversal plans, and diagnostics for efficient computation in the context of
𝓗² methods.

H2Trees builds upon ideas from [ClusterTrees](https://github.com/krcools/ClusterTrees.jl) and
provides several tree families:

* `TwoNTree`: an axis-aligned 2ⁿ-tree for organizing points or supported basis functions.
* `BlockTree`: a pair of test/trial trees for Petrov-Galerkin workflows.
* `BoundingBallTree`: a generic bounding-ball cluster tree.
* `KMeansTree`: a bounding-ball tree built from k-means clusters.
* `MetisTree` and `MetisForest`: graph-partitioned trees and forests for connected or disconnected geometries.
* `SimpleHybridTree`: a tree wrapper that separates upper and lower levels.

Tree construction uses a builder workflow:

```julia
using H2Trees

tree = buildtree(points; builder=TwoNTreeBuilder())
blocktree = buildtree(testpoints, trialpoints; builder=BlockTreeBuilder())
```

High-level builders provide sensible defaults. For example, `TwoNTreeBuilder()` builds an
adaptive tree with `minhalfsize=0` and `minvalues=70`, while BEAST spaces automatically resolve
to a protrusion-aware `TwoNTree` policy.

The package also provides:

* **Aggregation and disaggregation**: Plans for implementing aggregation and disaggregation algorithms for efficient computation.
* **Computation of translations**: Algorithms for computing translations between different tree levels.
* **Near/far diagnostics**: Tools such as `checkadmissibility` for validating assembled plans.
* **Plotting**: An interface to [PlotlyJS.jl](https://github.com/JuliaPlots/PlotlyJS.jl) for visualizing tree data structures.
* **Interface to BEAST**: An interface to the [BEAST](https://github.com/krcools/BEAST.jl) package for clustering of basis functions.
* **Interfaces to METIS and ParallelKMeans**: Optional extensions for graph- and clustering-based trees.

## Documentation

* Documentation for the [latest stable version](https://djukic14.github.io/H2Trees.jl/stable/).
* Documentation for the [development version](https://djukic14.github.io/H2Trees.jl/dev/).
