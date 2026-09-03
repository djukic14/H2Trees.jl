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

# A symmetry-reduced trait (storing one direction per orbit of the cubic lattice symmetries and
# recovering the rest through a stored symmetry ID) is the natural next member of this family.
# The geometry it needs already exists in `translations/latticesymmetry.jl` (`symmetrygroup`,
# `canonicalizetranslation`); what remains is threading a symmetry ID through `translationinfos`
# and into the consumer.
#
# Such a trait must be parameterized by the symmetry policy, and that policy must be the subgroup
# the CONSUMER can act with rather than the largest one the lattice admits: a representative
# found by minimizing over a larger group is in general unreachable from the actual offset by any
# element of the smaller one. See `canonicalizetranslation` for the reachability contract.

"""
    SymmetryDirectionInvariancePerLevel(policy) <: AbstractTranslationTrait

Share translation directions related by a lattice SYMMETRY, within each level.

Where [`DirectionInvariancePerLevel`](@ref) shares a direction only between pairs whose integer
offsets are equal, this shares one between every pair whose offsets lie in the same orbit of
`policy`. The stored direction is the orbit's canonical representative and each interaction
additionally records a symmetry ID, so the consumer can recover its own direction; see
[`canonicalizetranslation`](@ref) for the reachability contract and the direction convention.

`policy` MUST be the subgroup the consumer can act with, not the largest one the lattice admits.
A representative found by minimizing over a larger group is in general unreachable from the
actual offset by any element of a smaller one, and the mismatch is silent on offsets that happen
to be aligned. The default is [`OppositeSymmetry`](@ref), which asks only that the consumer's
representation be closed under the antipodal map.

Interactions carry a fourth field, `symmetryID`, that the other traits do not produce.
"""
struct SymmetryDirectionInvariancePerLevel{P<:AbstractLatticeSymmetryGroup} <:
       AbstractTranslationTrait
    policy::P
end

function SymmetryDirectionInvariancePerLevel()
    return SymmetryDirectionInvariancePerLevel(OppositeSymmetry())
end

"""
    DirectionInvariance <: AbstractTranslationTrait

Share equal translation directions across all relevant levels.

For generic trees, equality is approximate. For `TwoNTree`s, directions are
keyed by integer offsets on the finest common grid, allowing reuse across
levels.
"""
struct DirectionInvariance <: AbstractTranslationTrait end
