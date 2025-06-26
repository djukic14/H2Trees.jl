"""
    abstract type AbstractTranslationTrait

Abstract type for translation traits.
"""
abstract type AbstractTranslationTrait end

"""
    struct AllTranslations <: AbstractTranslationTrait

Represents the translation trait where all translations that occur are unique and are stored
individually.
"""
struct AllTranslations <: AbstractTranslationTrait end

#TODO: for very short translations the DirectionInvariancePerLevel exhibits numerical instability
# some translations might be counted twice or translations in the DisaggregationPlan are wrong
# and should not be performed: this is a bug

# translations on same level which have same direction are stored using same translation
"""
    struct DirectionInvariancePerLevel <: AbstractTranslationTrait

Represents the translation trait where translations on the same level with the same length
and direction are identical.
Therefore, only one version of the translation is stored.
"""
struct DirectionInvariancePerLevel <: AbstractTranslationTrait end

# TODO: for the MLFMA, for example, it is enough to consider translations in one quarter space and use symmetry
# at x=0,y=0,z=0 to get the other translations. (This is a case of "eierlegende Wollmilchsau".)
# struct QuarterSpaceSymmetryDirectionInvariancePerLevel <: AbstractTranslationTrait end

"""
    struct DirectionInvariance <: AbstractTranslationTrait

Represents the translation trait where the direction of translation is invariant,
i.e., translations in different directions are identical.
Therefore, only one version of the translation is stored.
"""
struct DirectionInvariance <: AbstractTranslationTrait end
