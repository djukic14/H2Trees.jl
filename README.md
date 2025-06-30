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

The [H2Trees](https://github.com/djukic14/H2Trees.jl) package provides a Julia implementation of tree data structures and algorithms for efficient computation in the context of 𝓗² methods.

H2Trees builds upon ideas from [ClusterTrees](https://github.com/krcools/ClusterTrees.jl) and provides a range of tree data structures, including:

* [`TwoNTree`](@ref): A 2ⁿ-tree data structure for organizing points in ℝⁿ.
* [`SimpleHybridTree`](@ref): A tree data structure that splits a tree
* [`BoundingBallTree`](@ref): A tree data structure using bounding balls.

The H2Trees package provides a range of features, including

* **Aggregation and disaggregation**: Plans for implementing aggregation and disaggregation algorithms for efficient computation.
* **Computation of translations**: Algorithms for computing translations between different tree levels.
* **Plotting**: An interface to [PlotlyJS.jl](https://github.com/JuliaPlots/PlotlyJS.jl) for visualizing tree data structures.
* **Interface to BEAST**: An interface to the [BEAST](https://github.com/krcools/BEAST.jl) package for clustering of basis functions.

## Documentation

* Documentation for the [latest stable version](https://djukic14.github.io/H2Trees.jl/stable/).
* Documentation for the [development version](https://djukic14.github.io/H2Trees.jl/dev/)
