# Translations

The following implementations of the [`AbstractTranslationTrait`](@ref) are available, enabling the specification of translations in a tree and reducing computational and storage requirements.

| Translation Trait                        | Description                                                                            |
|:-----------------------------------------|:---------------------------------------------------------------------------------------|
| [`AllTranslations`](@ref)                | Store and compute all translations individually.                                       |
| [`DirectionInvariance`](@ref)            | Treat translations with the same length and direction as identical.                    |
| [`DirectionInvariancePerLevel`](@ref)    | Treat translations on the same level with the same length and direction as identical.  |
| [`SymmetryDirectionInvariancePerLevel`](@ref) | Treat translations on the same level that are related by a lattice symmetry as identical. |

The translations can be computed with the [`translations`](@ref) function, which returns a plain
3-tuple `(translationinfos, translationdirections, relevantlevels)`:

- `translationinfos`: a vector of vectors containing `NamedTuple`s with fields `receivingnode`, `translatingnode`, and `translationID`, grouped by receiving level.
    The `translationID` is the id of the translation in the translation directions.
- `translationdirections`: the unique translation directions.
- `relevantlevels`: the levels covered by the translations.

## Symmetry-reduced translations

[`SymmetryDirectionInvariancePerLevel`](@ref) goes one step beyond
[`DirectionInvariancePerLevel`](@ref): rather than sharing a stored direction only between *exactly
equal* displacements, it shares one direction across a whole orbit of the lattice symmetry group, so
a single stored direction serves every displacement reachable from it by a symmetry. Which group is
used is set by the trait's `policy`: [`OppositeSymmetry`](@ref) (the conservative default, relating
`q` and `-q`), [`AxisPreservingSymmetry`](@ref), or [`FullLatticeSymmetry`](@ref). On the sphere
below this takes 436 stored directions to 218 and 24 respectively.

Its `translationinfos` records carry a **fourth field**, `symmetryID`, naming which group element
maps the stored (canonical) direction onto the pair's own. Every other trait produces three-field
records, so a consumer reading them positionally must know which trait produced them. The pair's own
direction is recovered as `applysymmetry(group[symmetryID], canonical)`.

The trait requires `TwoNTree` geometry, and refuses anything else in **two distinct ways** that
consumers should treat differently:

  - Any other tree, or a `BlockTree` with a non-`TwoNTree` side, raises an `ArgumentError`. A
    lattice symmetry has no meaning without a lattice, so no configuration of that geometry would
    work and nothing should catch this.
  - A `TwoNTree` pair whose two roots are offset by a vector that is not *itself* a lattice vector
    raises a catchable [`NonLatticeTranslationError`](@ref). Here the geometry is right and only this
    particular offset is not, so consumers are expected to catch it and fall back to
    [`DirectionInvariancePerLevel`](@ref), which deduplicates by displacement alone and is correct
    for any offset. The Petrov example below has coincident roots, so it takes the reduction.

### What the reduction looks like

The left panel is every translation one box has to receive from: its interaction list, the children
of its parent's neighbours minus its own neighbours. The right panel is what a symmetry-reduced
collection actually stores: one direction per orbit, each serving every arrow of its colour on the
left. The faint boxes on the right are the same interaction list, drawn for comparison.

```@raw html
<object data="../assets/plots/translation_symmetry_2d.html" type="text/html" style="width:100%; height:45vh;"> </object>
```

In three dimensions the same box has 189 far neighbours and the full lattice group has 48 elements,
which leaves 16 stored directions. The figure is interactive, so drag to rotate; the 189 offsets form
a hollow shell (a 7×7×7 block with the 3×3×3 near field removed), so a fixed viewpoint hides most of
it behind itself.

```@raw html
<object data="../assets/plots/translation_symmetry_3d.html" type="text/html" style="width:100%; height:55vh;"> </object>
```

Both figures use a box at child slot `(0,0)`/`(0,0,0)` and [`FullLatticeSymmetry`](@ref), and those
two choices are what make them exact rather than indicative. Canonical representatives are chosen
globally, so for another child slot some of them would not be offsets of the box being drawn and the
right panel would contain arrows the left never showed. And under the full group one box's list
already meets every orbit the whole level has, so 7 and 16 are simultaneously this box's counts and
the level-wide stored counts. Under [`OppositeSymmetry`](@ref) those differ (19 against 20 in two
dimensions), and the picture would need an asterisk.

This can for example look like

```@example translations
using CompScienceMeshes # hide
using H2Trees # hide

m = meshsphere(1.0, 0.1)
tree = buildtree(vertices(m); builder=TwoNTreeBuilder(; minhalfsize=0.1))

# we are going to aggregate all nodes, even though that might not be needed, and we use the
# AggregateMode; the default isnear() function decides which translations occur
plans = H2Trees.buildplans(
    tree; builder=H2Trees.PlanBuilder(; aggregatenode=H2Trees.AggregateAllNodesFunctor())
)

# the translations can be found in the testdisaggregationplan in AggregateMode
for translationtrait in [H2Trees.AllTranslations(), H2Trees.DirectionInvariance(),  H2Trees.DirectionInvariancePerLevel(),
                         H2Trees.SymmetryDirectionInvariancePerLevel(H2Trees.OppositeSymmetry()),
                         H2Trees.SymmetryDirectionInvariancePerLevel(H2Trees.FullLatticeSymmetry())]
    infos, directions, relevantlevels = H2Trees.translations(
        tree, plans.testdisaggregationplan, translationtrait
    )
    println("For translationtrait $(typeof(translationtrait)) we have $(length(directions)) unique translations.")
end
```

and in the Petrov-Galerkin case

```@example translations2
using CompScienceMeshes # hide
using H2Trees # hide

mx = meshsphere(1.0, 0.1)
my = meshsphere(2.0, 0.1)

tree = buildtree(vertices(mx), vertices(my); builder=BlockTreeBuilder(; test=TwoNTreeBuilder(; minhalfsize=0.1), trial=TwoNTreeBuilder(; minhalfsize=0.1)))
testtree = H2Trees.testtree(tree)
trialtree = H2Trees.trialtree(tree)

# we are going to aggregate all nodes, even though that might not be needed, and we use the
# AggregateMode; the default isnear() function decides which translations occur; buildplans
# dispatches to the Petrov path automatically since `tree` is a BlockTree
plans = H2Trees.buildplans(
    tree; builder=H2Trees.PlanBuilder(; aggregatenode=H2Trees.AggregateAllNodesFunctor())
)

# the translations can be found in the testdisaggregationplan in AggregateMode
for translationtrait in [H2Trees.AllTranslations(), H2Trees.DirectionInvariance(),  H2Trees.DirectionInvariancePerLevel(),
                         H2Trees.SymmetryDirectionInvariancePerLevel(H2Trees.OppositeSymmetry()),
                         H2Trees.SymmetryDirectionInvariancePerLevel(H2Trees.FullLatticeSymmetry())]
    infos, directions, relevantlevels = H2Trees.translations(
        tree, plans.testdisaggregationplan, translationtrait
    )
    println("For translationtrait $(typeof(translationtrait)) we have $(length(directions)) unique translations.")
end
```
