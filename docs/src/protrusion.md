# Protrusion Policy

Protrusion measures how far a stored value extends beyond the box it is assigned to, normalized
by `2 * halfsize`. It only matters for [`TwoNTree`](@ref): points are dimensionless, but BEAST
basis functions have physical support that can straddle a box boundary. `TwoNTreeBuilder`'s
`protrusion` option controls whether subdivision stops early to keep that support from spilling
too far out of its box.

## Policies

- [`NoProtrusionCheck`](@ref) — never stop early for protrusion. Default for plain point
  collections, which have no extent.
- [`AutoProtrusionCheck`](@ref) — let the input type choose (the default). Plain points resolve
  to `NoProtrusionCheck()`; BEAST spaces resolve to `ProtrusionCheck(; max=0.25, compute=BEASTProtrusionFunctor(space))`.
- [`ProtrusionCheck`](@ref)`(; max, compute=ComputeProtrusionFunctor())` — reject a split when any
  value protrudes at least `max` (normalized). `compute(center, halfsize, value)` returns the
  normalized protrusion; `max=0.25` means a value may extend at most `0.5 * halfsize` beyond the box.

## Overriding

```@example protrusion1
using BEAST, CompScienceMeshes # hide
using H2Trees

m = meshsphere(1.0, 0.1)
X = raviartthomas(m)

default = TwoNTree(X)                                           # auto -> max=0.25
looser = TwoNTree(X; builder=TwoNTreeBuilder(; protrusion=ProtrusionCheck(; max=0.5)))
off = TwoNTree(X; builder=TwoNTreeBuilder(; protrusion=NoProtrusionCheck()))
nothing #hide
```

A custom `compute` functor just needs to be callable as `(center, halfsize, value) -> Float64`.

## Diagnosing an existing tree

[`H2Trees.protrusionreport`](@ref) reports the worst level; [`H2Trees.levelprotrusions`](@ref)
returns the full per-level vector. A level at or above `0.5` means values are protruding by a
full box or more, and a warning is emitted automatically:

```@example protrusion1
H2Trees.protrusionreport(default)
```
