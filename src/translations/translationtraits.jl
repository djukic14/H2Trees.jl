abstract type AbstractTranslationTrait end

struct AllTranslations <: AbstractTranslationTrait end

#TODO: for very short translations the DirectionInvariancePerLevel exhibits numerical instability
# some translations might be counted twice or translations in the DisaggregationPlan are wrong
# and should not be performed: this is a bug

# translations on same level which have same direction are stored using same translation
struct DirectionInvariancePerLevel <: AbstractTranslationTrait end

# TODO: for the MLFMA, for example, it is enough to consider translations in one quarter space and use symmetry
# at x=0,y=0,z=0 to get the other translations. (This is a case of "eierlegende Wollmilchsau".)
# struct QuarterSpaceSymmetryDirectionInvariancePerLevel <: AbstractTranslationTrait end

struct DirectionInvariance <: AbstractTranslationTrait end
