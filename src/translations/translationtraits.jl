"""
    AbstractTranslationTrait

Supertype for translation-direction storage traits.

Translation traits control how [`translations`](@ref) deduplicates scheduled
translation vectors before downstream translation operators are built.
"""
abstract type AbstractTranslationTrait end

"""
    AllTranslations <: AbstractTranslationTrait

Store one translation direction for every scheduled pair.

Use this when no direction reuse is available or when exact pair identity matters.
"""
struct AllTranslations <: AbstractTranslationTrait end

"""
    DirectionInvariancePerLevel <: AbstractTranslationTrait

Share equal translation directions within each level.

For generic trees, directions are deduplicated with approximate comparison. For
`TwoNTree`s, the specialized translation path uses integer grid offsets, so
directions on the same level are keyed exactly.
"""
struct DirectionInvariancePerLevel <: AbstractTranslationTrait end

# A future MLFMA-oriented trait could store only a symmetry-reduced subset of
# directions and recover the remaining ones by coordinate symmetries.
# struct QuarterSpaceSymmetryDirectionInvariancePerLevel <: AbstractTranslationTrait end

"""
    DirectionInvariance <: AbstractTranslationTrait

Share equal translation directions across all relevant levels.

For generic trees, equality is approximate. For `TwoNTree`s, directions are
keyed by integer offsets on the finest common grid, allowing reuse across
levels.
"""
struct DirectionInvariance <: AbstractTranslationTrait end
