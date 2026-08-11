"""
    SEBB

Smallest Enclosing Ball of Balls (SEBB) in dimensions 1, 2, and 3.

Exact (the certified minimal ball) whenever the configuration is resolvable at the caller's
precision. A degenerate input never raises: the public solver warns once and returns a strictly
enclosing approximation instead. Tests and diagnostics can inspect the internal certificate
helper to distinguish exact from fallback results.

This submodule is tree-independent. H2Trees adapts trees to SEBB, not the reverse.

# Public API

  - [`Ball`](@ref), [`center`](@ref), [`radius`](@ref)
  - [`encloses`](@ref)
  - [`smallest_enclosing_ball`](@ref)
"""
module SEBB

using LinearAlgebra
using StaticArrays

export Ball
export center, radius
export encloses
export smallest_enclosing_ball

include("ball.jl")
include("tolerances.jl")
include("predicates.jl")
include("roots.jl")
include("support.jl")
include("fallback.jl")
include("enumeration.jl")
include("validation.jl")

end # module SEBB
