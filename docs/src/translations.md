# Translations

The following implementations of the [`AbstractTranslationTrait`](@ref) are available, enabling the specification of translations in a tree and reducing computational and storage requirements.

| Translation Trait                        | Description                                                                            |
|:-----------------------------------------|:---------------------------------------------------------------------------------------|
| [`AllTranslations`](@ref)                | Store and compute all translations individually.                                       |
| [`DirectionInvariance`](@ref)            | Treat translations with the same length and direction as identical.                    |
| [`DirectionInvariancePerLevel`](@ref)    | Treat translations on the same level with the same length and direction as identical.  |

The translations can be computed with the [`translations`](@ref) function, which returns a plain
3-tuple `(translationinfos, translationdirections, relevantlevels)`:

- `translationinfos`: a vector of vectors containing `NamedTuple`s with fields `receivingnode`, `translatingnode`, and `translationID`, grouped by receiving level.
    The `translationID` is the id of the translation in the translation directions.
- `translationdirections`: the unique translation directions.
- `relevantlevels`: the levels covered by the translations.

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
for translationtrait in [H2Trees.AllTranslations(), H2Trees.DirectionInvariance(),  H2Trees.DirectionInvariancePerLevel()]
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
for translationtrait in [H2Trees.AllTranslations(), H2Trees.DirectionInvariance(),  H2Trees.DirectionInvariancePerLevel()]
    infos, directions, relevantlevels = H2Trees.translations(
        tree, plans.testdisaggregationplan, translationtrait
    )
    println("For translationtrait $(typeof(translationtrait)) we have $(length(directions)) unique translations.")
end
```
