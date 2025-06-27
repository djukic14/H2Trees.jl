# Translations

The following implementations of the [`AbstractTranslationTrait`](@ref) are available, enabling the specification of translations in a tree and reducing computational and storage requirements.

| Translation Trait                        | Description                                                                            |
|:-----------------------------------------|:---------------------------------------------------------------------------------------|
| [`AllTranslations`](@ref)                | Store and compute all translations individually.                                       |
| [`DirectionInvariance`](@ref)            | Treat translations with the same length and direction as identical.                    |
| [`DirectionInvariancePerLevel`](@ref)    | Treat translations on the same level with the same length and direction as identical.  |

The translations can be computed with the [`translations`](@ref) function.
The result is a tuple containing two vectors:

- The first vector contains `NamedTuple`s with fields `receivingnode`, `translatingnode`, and `translationID`.
    The `translationID` is the id of the translation in the translation directions.
- The second vector contains the translation directions.

This can for example look like

```@example translations
using CompScienceMeshes # hide
using H2Trees # hide

m = meshsphere(1.0, 0.1)
tree = TwoNTree(vertices(m), 0.1)

# define iterator, which specifies which translations occur (we use the default isnear() function) 
tfiterator = H2Trees.TranslatingNodesIterator(; isnear=H2Trees.isnear())
# we are going to aggregate all nodes, even though that might not be needed and we use the AggregateMode
plans = H2Trees.galerkinplans(tree, H2Trees.AggregateAllNodesFunctor(), tfiterator, H2Trees.AggregateMode())

# the translations can be found in the testdisaggregationplan in AggregateMode
for translationtrait in [H2Trees.AllTranslations(), H2Trees.DirectionInvariance(),  H2Trees.DirectionInvariancePerLevel()]
    translationinfo, translations = H2Trees.translations(tree, plans.testdisaggregationplan, translationtrait)
    println("For translationtrait $(typeof(translationtrait)) we have $(length(translations)) unique translations.")
end
```

and in the Petrov-Galerkin case

```@example translations2
using CompScienceMeshes # hide
using H2Trees # hide

mx = meshsphere(1.0, 0.1)
my = meshsphere(2.0, 0.1)

tree = TwoNTree(vertices(mx), vertices(my), 0.1)
testtree = H2Trees.testtree(tree)
trialtree = H2Trees.trialtree(tree)

# define iterator, which specifies which translations occur (we use the default isnear() function) 
tfiterator = H2Trees.TranslatingNodesIterator(; isnear=H2Trees.isnear())

# we are going to aggregate all nodes, even though that might not be needed and we use the AggregateMode
plans = H2Trees.petrovplans(tree, H2Trees.AggregateAllNodesFunctor(), tfiterator, H2Trees.AggregateMode())

# the translations can be found in the testdisaggregationplan in AggregateMode
for translationtrait in [H2Trees.AllTranslations(), H2Trees.DirectionInvariance(),  H2Trees.DirectionInvariancePerLevel()]
    translationinfo, translations = H2Trees.translations(tree, plans.testdisaggregationplan, translationtrait)
    println("For translationtrait $(typeof(translationtrait)) we have $(length(translations)) unique translations.")
end
```
