# 𝓗² Trees

The [H2Trees](https://github.com/djukic14/H2Trees.jl) package provides a Julia implementation of tree data structures and algorithms for efficient computation in the context of 𝓗² methods.

H2Trees builds upon ideas from [ClusterTrees](https://github.com/krcools/ClusterTrees.jl) and provides a range of tree data structures, including:

* [`TwoNTree`](@ref): A 2ⁿ-tree data structure for organizing points in ℝⁿ.
* [`SimpleHybridTree`](@ref): A tree data structure that splits a tree into two parts for hybrid methods.
* [`BoundingBallTree`](@ref): A tree data structure using bounding balls.

See [Tree Families](tree_families.md) for a full comparison and how to pick one, and
[Builder Workflow](builders.md) for the construction API shared by all of them.

The H2Trees package provides a range of features, including

* **Aggregation and disaggregation**: Plans for implementing aggregation and disaggregation algorithms for efficient computation. See [Plans](plans/plans.md) for details.
* **Computation of translations**: Algorithms for computing translations between different tree levels.
* **Plotting**: An interface to [PlotlyJS.jl](https://github.com/JuliaPlots/PlotlyJS.jl) for visualizing tree data structures.
* **Interface to BEAST**: An interface to the [BEAST](https://github.com/krcools/BEAST.jl) package for clustering of basis functions.
* **Interface to METIS**: An interface to [Metis.jl](https://github.com/JuliaSparse/Metis.jl) for graph-based partitioning strategies.
* **Interface to ParallelKMeans**: An interface to [ParallelKMeans.jl](https://github.com/PyDataBlog/ParallelKMeans.jl) for k-means based partitioning strategies.

## Goals

The primary goal of the [H2Trees](https://github.com/djukic14/H2Trees.jl) package is to provide an efficient and flexible framework for computing trees for 𝓗² methods, and related quantities.
The package aims to provide a range of features and tools for working with tree data structures, including aggregation, disaggregation, and plotting.

## Related Packages

The [H2Trees](https://github.com/djukic14/H2Trees.jl) package is related to the following packages:

* [ClusterTrees](https://github.com/krcools/ClusterTrees.jl)
* [BEAST](https://github.com/krcools/BEAST.jl)
* [H2Factory](https://github.com/djukic14/H2Factory.jl)
* [MLFMA](https://github.com/djukic14/MLFMA.jl)
* [iFMM](https://github.com/djukic14/iFMM.jl)
* [HybridFMM](https://github.com/djukic14/HybridFMM.jl)

## Documentation

This documentation provides an overview of the [H2Trees](https://github.com/djukic14/H2Trees.jl) package, including its features, goals, and related packages. For more information, please see the individual documentation pages for each module and function.
