"""
    _CenterFunctor(tree)

Callable wrapper returning `center(tree, node)`.
"""
struct _CenterFunctor{T}
    tree::T
end

function (f::_CenterFunctor)(node::Int)
    return center(f.tree, node)
end

"""
    _HalfSizeFunctor(tree)

Callable wrapper returning `halfsize(tree, node)`.
"""
struct _HalfSizeFunctor{T}
    tree::T
end

function (f::_HalfSizeFunctor)(node::Int)
    return halfsize(f.tree, node)
end

"""
    _LevelFunctor(tree)

Callable wrapper returning `level(tree, node)`.
"""
struct _LevelFunctor{T}
    tree::T
end

function (f::_LevelFunctor)(node::Int)
    return level(f.tree, node)
end

"""
    foreachtranslationpair(f, translatingplan, relevantlevels, receivinglevel)

Call `f(levelid, receivingnode, translatingnode)` for every pair scheduled by a
translating plan.

`levelid` indexes `relevantlevels` according to the receiving node's level.
"""
function foreachtranslationpair(
    f, translatingplan::AbstractPlan, relevantlevels, receivinglevel
)
    relevantlevelsdict = Dict(zip(relevantlevels, collect(eachindex(relevantlevels))))

    for level in relevantlevels
        for receivingnode in receivingnodes(translatingplan, level)
            rlevelid = relevantlevelsdict[receivinglevel(receivingnode)]
            for translatingnode in translatingplan[receivingnode, level]
                f(rlevelid, receivingnode, translatingnode)
            end
        end
    end

    return nothing
end

"""
    _translationinfo(receivingnode, translatingnode, translationID)

Create the compact metadata record stored in `translationinfos`.
"""
function _translationinfo(receivingnode::Int, translatingnode::Int, translationID::Int)
    return (
        receivingnode=receivingnode,
        translatingnode=translatingnode,
        translationID=translationID,
    )
end

"""
    settranslationinfo!(translationinfos, rlevelid, levelindex, receivingnode,
        translatingnode, translationID)

Store one translation metadata record at a known level-local index.
"""
function settranslationinfo!(
    translationinfos,
    rlevelid::Int,
    levelindex::Int,
    receivingnode::Int,
    translatingnode::Int,
    translationID::Int,
)
    translationinfos[rlevelid][levelindex] = _translationinfo(
        receivingnode, translatingnode, translationID
    )
    return nothing
end

"""
    appendtranslationinfo!(translationinfos, rlevelid, receivingnode,
        translatingnode, translationID)

Append one translation metadata record to a level-local vector.
"""
function appendtranslationinfo!(
    translationinfos,
    rlevelid::Int,
    receivingnode::Int,
    translatingnode::Int,
    translationID::Int,
)
    push!(
        translationinfos[rlevelid],
        _translationinfo(receivingnode, translatingnode, translationID),
    )
    return nothing
end

"""
    _symmetrytranslationinfo(receivingnode, translatingnode, translationID, symmetryID)

Create the metadata record stored by [`SymmetryDirectionInvariancePerLevel`](@ref).

It carries a fourth field the other traits do not produce, rather than widening the common
record: every existing consumer reads a three-field record, and adding a field they would ignore
to all of them would change a type used across package boundaries for the benefit of one trait.
"""
function _symmetrytranslationinfo(
    receivingnode::Int, translatingnode::Int, translationID::Int, symmetryID::Int
)
    return (
        receivingnode=receivingnode,
        translatingnode=translatingnode,
        translationID=translationID,
        symmetryID=symmetryID,
    )
end

"""
    appendsymmetrytranslationinfo!(translationinfos, rlevelid, receivingnode,
        translatingnode, translationID, symmetryID)

Append one symmetry-aware translation metadata record to a level-local vector.
"""
function appendsymmetrytranslationinfo!(
    translationinfos,
    rlevelid::Int,
    receivingnode::Int,
    translatingnode::Int,
    translationID::Int,
    symmetryID::Int,
)
    push!(
        translationinfos[rlevelid],
        _symmetrytranslationinfo(receivingnode, translatingnode, translationID, symmetryID),
    )
    return nothing
end

"""
    translations(tree, translatingplan, translationtrait)

Compute translation vectors for a translating plan.

Returns `(translationinfos, translationdirections, relevantlevels)`.
`translationinfos[levelid]` stores named tuples with `receivingnode`,
`translatingnode`, and `translationID`; `translationdirections[translationID]`
stores the vector used by that pair.

The `translationtrait` controls how directions are deduplicated:

  - [`AllTranslations`](@ref): every scheduled pair gets its own direction.

  - [`DirectionInvariancePerLevel`](@ref): equal directions are shared within a level.

  - [`DirectionInvariance`](@ref): equal directions are shared across all relevant levels.

  - [`SymmetryDirectionInvariancePerLevel`](@ref): directions are shared within a level up to a
    LATTICE SYMMETRY, so one stored direction serves a whole orbit rather than only exact
    duplicates. Its `translationinfos` carry an extra `symmetryID` naming which group element maps
    the stored (canonical) direction onto the pair's own: every other trait's records have three
    fields, these have four, so a consumer reading them positionally must know which trait produced
    them.

    TWO REFUSALS ARE PART OF ITS CONTRACT, and they are deliberately different types.
    `TwoNTree` geometry is REQUIRED: any other tree, or a `BlockTree` with a non-`TwoNTree` side, is
    an `ArgumentError`, because a lattice symmetry has no meaning without a lattice. A `TwoNTree`
    pair whose roots are offset by a non-lattice vector instead throws a catchable
    [`NonLatticeTranslationError`](@ref), which consumers are expected to handle by falling back to
    `DirectionInvariancePerLevel`. Configuration error versus fallback condition.

For `BlockTree`s, receiving and translating sides are resolved from the
translating plan.
"""
function translations(tree, translatingplan::AbstractPlan, translationtrait)
    if !istranslatingplan(translatingplan)
        throw(ArgumentError("translations require a translating plan"))
    end
    return translations(tree, treetrait(tree), translatingplan, translationtrait)
end

function translations(
    tree, ::AbstractTreeTrait, translatingplan::AbstractPlan, translationtrait
)
    relevantlevels = mintranslationlevel(translatingplan):(H2Trees.levels(tree)[end])

    return (
        _translations(
            translatingplan,
            relevantlevels,
            _CenterFunctor(tree),
            _CenterFunctor(tree),
            _LevelFunctor(tree),
            Val{eltype(tree)}(),
            translationtrait,
        )...,
        relevantlevels,
    )
end

"""
Refuse [`SymmetryDirectionInvariancePerLevel`](@ref) on a geometry that cannot carry it.

AN EXPLICIT CONTRACT, NOT AN ACCIDENTAL DISPATCH GAP. Without this method the generic
`translations(tree, ::AbstractTreeTrait, ...)` above forwards to a seven-argument `_translations`,
while the symmetry trait's own method takes eight (it additionally needs `receivinghalfsize`). The
result was a `MethodError` naming `_translations` and a list of functors, which says nothing about
the trait, the tree, or what to use instead.

DISTINCT FROM [`NonLatticeTranslationError`](@ref), and the two must not be conflated. That one is a
CATCHABLE FALLBACK CONDITION: the geometry is right, this particular root offset is not, and a
consumer is expected to catch it and drop to `DirectionInvariancePerLevel`. This is a programmer or
configuration error (the geometry itself cannot support lattice symmetries at any offset), so it
is an `ArgumentError` and nothing should be catching it.

The `::isTwoNTree` method below is strictly more specific and still wins.
"""
function translations(
    tree,
    ::AbstractTreeTrait,
    translatingplan::AbstractPlan,
    ::SymmetryDirectionInvariancePerLevel,
)
    return throw(
        ArgumentError(
            "SymmetryDirectionInvariancePerLevel requires TwoNTree lattice geometry, but this " *
            "tree is a $(typeof(tree)) ($(treetrait(tree))). A lattice symmetry maps a " *
            "displacement to another displacement on the SAME lattice, which needs node centers " *
            "separated by integer multiples of the level halfsize: a property of the regular " *
            "2^D subdivision, not of trees in general. Use `DirectionInvariancePerLevel`, which " *
            "deduplicates by displacement alone and is correct for any geometry.",
        ),
    )
end

function translations(
    tree, ::isTwoNTree, translatingplan::AbstractPlan, ::DirectionInvariancePerLevel
)
    relevantlevels = mintranslationlevel(translatingplan):(levels(tree)[end])

    return (
        _translations(
            translatingplan,
            relevantlevels,
            _CenterFunctor(tree),
            _CenterFunctor(tree),
            _HalfSizeFunctor(tree),
            _LevelFunctor(tree),
            Val{eltype(tree)}(),
            DirectionInvariancePerLevel(),
        )...,
        relevantlevels,
    )
end

function translations(
    tree, ::isTwoNTree, translatingplan::AbstractPlan, ::DirectionInvariance
)
    relevantlevels = mintranslationlevel(translatingplan):(levels(tree)[end])

    return (
        _translations(
            translatingplan,
            relevantlevels,
            _CenterFunctor(tree),
            _CenterFunctor(tree),
            minhalfsize(tree),
            _LevelFunctor(tree),
            Val{eltype(tree)}(),
            DirectionInvariance(),
        )...,
        relevantlevels,
    )
end

function translations(
    tree,
    ::isTwoNTree,
    translatingplan::AbstractPlan,
    translationtrait::SymmetryDirectionInvariancePerLevel,
)
    relevantlevels = mintranslationlevel(translatingplan):(levels(tree)[end])

    return (
        _translations(
            translatingplan,
            relevantlevels,
            _CenterFunctor(tree),
            _CenterFunctor(tree),
            _HalfSizeFunctor(tree),
            _LevelFunctor(tree),
            Val{eltype(tree)}(),
            translationtrait,
        )...,
        relevantlevels,
    )
end

function translations(tree, ::isBlockTree, translatingplan::AbstractPlan, translationtrait)
    receive = receivingtree(tree, translatingplan)
    translate = translatingtree(tree, translatingplan)
    return translations(
        receive,
        translate,
        treetrait(receive),
        treetrait(translate),
        translatingplan,
        translationtrait,
    )
end

# THE SAME FORWARDING, RESTATED FOR THE SYMMETRY TRAIT ALONE, and it has to exist rather than being
# covered by the method above. Against the geometry refusal
# `(::Any, ::AbstractTreeTrait, ::AbstractPlan, ::SymmetryDirectionInvariancePerLevel)` the method
# above is more specific in argument 2 and less specific in argument 4, so neither wins: a
# `BlockTree` asking for a symmetry is AMBIGUOUS between forwarding and refusing. Julia resolved it
# by refusing, which broke every block-tree symmetry case, including the ones that must work.
#
# This method is strictly more specific than both, so it settles the ambiguity in favour of
# forwarding. The two sides' own traits are then what decide, in the six-argument methods below:
# `TwoNTree`/`TwoNTree` is served, anything else is refused there. A `BlockTree` is not itself
# lattice or non-lattice: only its sides are, which is precisely why the decision belongs one
# level down and not here.
function translations(
    tree,
    ::isBlockTree,
    translatingplan::AbstractPlan,
    translationtrait::SymmetryDirectionInvariancePerLevel,
)
    receive = receivingtree(tree, translatingplan)
    translate = translatingtree(tree, translatingplan)
    return translations(
        receive,
        translate,
        treetrait(receive),
        treetrait(translate),
        translatingplan,
        translationtrait,
    )
end

function translations(
    receivetree,
    translatingtree,
    ::AbstractTreeTrait,
    ::AbstractTreeTrait,
    translatingplan::AbstractPlan,
    translationtrait,
)
    leaflevel = min(levels(receivetree)[end], levels(translatingtree)[end])
    relevantlevels = mintranslationlevel(translatingplan):leaflevel

    return (
        _translations(
            translatingplan,
            relevantlevels,
            _CenterFunctor(receivetree),
            _CenterFunctor(translatingtree),
            _LevelFunctor(receivetree),
            Val{promote_type(eltype(receivetree), eltype(translatingtree))}(),
            translationtrait,
        )...,
        relevantlevels,
    )
end

"""
The block-tree twin of the single-tree refusal above: BOTH sides must be `TwoNTree`s.

One side being non-lattice is enough to break the argument: the symmetry maps a displacement
between the two trees, so it needs both endpoints on a common lattice. The `::isTwoNTree, ::isTwoNTree` method below is strictly more specific and still wins; this catches every mixed pair.
"""
function translations(
    receivetree,
    translatingtree,
    ::AbstractTreeTrait,
    ::AbstractTreeTrait,
    translatingplan::AbstractPlan,
    ::SymmetryDirectionInvariancePerLevel,
)
    return throw(
        ArgumentError(
            "SymmetryDirectionInvariancePerLevel requires BOTH block-tree sides to be TwoNTrees, " *
            "but the receiving side is a $(typeof(receivetree)) ($(treetrait(receivetree))) and " *
            "the translating side a $(typeof(translatingtree)) " *
            "($(treetrait(translatingtree))). The symmetry maps a displacement BETWEEN the two " *
            "trees, so it needs both endpoints on one common lattice; one non-lattice side is " *
            "enough to break that. Use `DirectionInvariancePerLevel`, which deduplicates by " *
            "displacement alone and is correct for any pair of geometries.",
        ),
    )
end

function translations(
    receivetree,
    translatingtree,
    ::isTwoNTree,
    ::isTwoNTree,
    translatingplan::AbstractPlan,
    ::DirectionInvariancePerLevel,
)
    leaflevel = min(levels(receivetree)[end], levels(translatingtree)[end])
    relevantlevels = mintranslationlevel(translatingplan):leaflevel

    return (
        _translations(
            translatingplan,
            relevantlevels,
            _CenterFunctor(receivetree),
            _CenterFunctor(translatingtree),
            _HalfSizeFunctor(receivetree),
            _LevelFunctor(receivetree),
            Val{promote_type(eltype(receivetree), eltype(translatingtree))}(),
            DirectionInvariancePerLevel();
            offset=H2Trees.center(receivetree, H2Trees.root(receivetree)) -
                   H2Trees.center(translatingtree, H2Trees.root(translatingtree)),
        )...,
        relevantlevels,
    )
end

function translations(
    receivetree,
    translatingtree,
    ::isTwoNTree,
    ::isTwoNTree,
    translatingplan::AbstractPlan,
    translationtrait::SymmetryDirectionInvariancePerLevel,
)
    leaflevel = min(levels(receivetree)[end], levels(translatingtree)[end])
    relevantlevels = mintranslationlevel(translatingplan):leaflevel

    return (
        _translations(
            translatingplan,
            relevantlevels,
            _CenterFunctor(receivetree),
            _CenterFunctor(translatingtree),
            _HalfSizeFunctor(receivetree),
            _LevelFunctor(receivetree),
            Val{promote_type(eltype(receivetree), eltype(translatingtree))}(),
            translationtrait;
            offset=H2Trees.center(receivetree, H2Trees.root(receivetree)) -
                   H2Trees.center(translatingtree, H2Trees.root(translatingtree)),
        )...,
        relevantlevels,
    )
end

function translations(
    receivetree,
    translatingtree,
    ::isTwoNTree,
    ::isTwoNTree,
    translatingplan::AbstractPlan,
    ::DirectionInvariance,
)
    leaflevel = min(levels(receivetree)[end], levels(translatingtree)[end])
    relevantlevels = mintranslationlevel(translatingplan):leaflevel

    return (
        _translations(
            translatingplan,
            relevantlevels,
            _CenterFunctor(receivetree),
            _CenterFunctor(translatingtree),
            min(minhalfsize(receivetree), minhalfsize(translatingtree)),
            _LevelFunctor(receivetree),
            Val{promote_type(eltype(receivetree), eltype(translatingtree))}(),
            DirectionInvariance();
            offset=center(receivetree, root(receivetree)) -
                   center(translatingtree, root(translatingtree)),
        )...,
        relevantlevels,
    )
end

"""
    _translations(..., ::AllTranslations)

Materialize one translation direction per scheduled pair.
"""
function _translations(
    translatingplan::AbstractPlan,
    relevantlevels,
    receivecenter,
    translatingcenter,
    receivinglevel,
    ::Val{ELTYPE},
    ::AllTranslations,
) where {ELTYPE}
    # count number of translations on each level
    ntranslationsperlevel = zeros(Int, length(relevantlevels))
    foreachtranslationpair(
        translatingplan, relevantlevels, receivinglevel
    ) do rlevelid, _, _
        ntranslationsperlevel[rlevelid] += 1
        return nothing
    end

    # allocate required memory
    translationinfos = Vector{
        Vector{@NamedTuple{receivingnode::Int,translatingnode::Int,translationID::Int}}
    }(
        undef, length(relevantlevels)
    )
    for levelid in eachindex(relevantlevels)
        translationinfos[levelid] = Vector{
            @NamedTuple{receivingnode::Int,translatingnode::Int,translationID::Int}
        }(
            undef, ntranslationsperlevel[levelid]
        )
    end
    translations = Vector{ELTYPE}(undef, sum(ntranslationsperlevel))

    # compute translations
    translationID = 1
    translationIDlevel = ones(Int, length(relevantlevels))
    foreachtranslationpair(
        translatingplan, relevantlevels, receivinglevel
    ) do rlevelid, receivingnode, translatingnode
        translation = receivecenter(receivingnode) - translatingcenter(translatingnode)

        translations[translationID] = translation

        settranslationinfo!(
            translationinfos,
            rlevelid,
            translationIDlevel[rlevelid],
            receivingnode,
            translatingnode,
            translationID,
        )

        translationID += 1
        translationIDlevel[rlevelid] += 1
        return nothing
    end

    return translationinfos, translations
end

"""
    _translations(..., ::DirectionInvariancePerLevel)

Deduplicate equal translation vectors separately on each relevant level.
"""
function _translations(
    translatingplan::AbstractPlan,
    relevantlevels,
    receivecenter,
    translatingcenter,
    receivinglevel,
    ::Val{ELTYPE},
    ::DirectionInvariancePerLevel;
    isapprox=Base.isapprox,
) where {ELTYPE}
    translations = [ELTYPE[] for _ in relevantlevels]
    translationIDs = [Int[] for _ in relevantlevels]
    translationinfos = [
        @NamedTuple{receivingnode::Int, translatingnode::Int, translationID::Int}[] for
        _ in relevantlevels
    ]

    translationID = 1
    temptranslationID = 0

    foreachtranslationpair(
        translatingplan, relevantlevels, receivinglevel
    ) do rlevelid, receivingnode, translatingnode
        translation = receivecenter(receivingnode) - translatingcenter(translatingnode)
        translationindexinarray = findfirst(isapprox(translation), translations[rlevelid])

        if isnothing(translationindexinarray)
            push!(translations[rlevelid], translation)
            push!(translationIDs[rlevelid], translationID)
            temptranslationID = translationID
            translationID += 1
        else
            temptranslationID = translationIDs[rlevelid][translationindexinarray]
        end

        appendtranslationinfo!(
            translationinfos, rlevelid, receivingnode, translatingnode, temptranslationID
        )
        return nothing
    end

    translationarray = Vector{ELTYPE}(undef, translationID - 1)

    for i in eachindex(translations)
        for j in eachindex(translations[i])
            translationarray[translationIDs[i][j]] = translations[i][j]
        end
    end

    return translationinfos, translationarray
end

"""
    _translations(..., ::DirectionInvariance)

Deduplicate equal translation vectors across all relevant levels.
"""
function _translations(
    translatingplan::AbstractPlan,
    relevantlevels,
    receivecenter,
    translatingcenter,
    receivinglevel,
    ::Val{ELTYPE},
    ::DirectionInvariance;
    isapprox=Base.isapprox,
) where {ELTYPE}
    translations = ELTYPE[]
    translationIDs = Int[]
    translationinfos = [
        @NamedTuple{receivingnode::Int, translatingnode::Int, translationID::Int}[] for
        _ in relevantlevels
    ]

    translationID = 1
    temptranslationID = 0

    foreachtranslationpair(
        translatingplan, relevantlevels, receivinglevel
    ) do rlevelid, receivingnode, translatingnode
        translation = receivecenter(receivingnode) - translatingcenter(translatingnode)
        translationindexinarray = findfirst(isapprox(translation), translations)

        if isnothing(translationindexinarray)
            push!(translations, translation)
            push!(translationIDs, translationID)
            temptranslationID = translationID
            translationID += 1
        else
            temptranslationID = translationIDs[translationindexinarray]
        end

        appendtranslationinfo!(
            translationinfos, rlevelid, receivingnode, translatingnode, temptranslationID
        )
        return nothing
    end
    return translationinfos, translations
end

# The displacement in units of the level halfsize, REQUIRED to be a lattice vector.
#
# THIS IS WHERE THE SYMMETRY PATH DIVERGES FROM THE PLAIN ONE, and it is not a stylistic choice.
# `DirectionInvariancePerLevel` subtracts the block tree's root-center `offset` before rounding,
# which is sound there: every interaction sharing a lattice index has the SAME physical direction
# `halfsize * index + offset`, so one stored translation serves them all.
#
# A symmetry merges DIFFERENT indices. The stored direction is then
# `halfsize * canonical + offset` while the interaction needs
# `halfsize * (S * canonical) + offset`, and applying `S` to the stored one gives
# `halfsize * (S * canonical) + S * offset`. Those agree only when `S * offset == offset`. For a
# generic offset that is the identity alone, so with a displaced pair of roots the reduction is
# not merely less effective, it is WRONG, and wrong silently: the shapes match and the field is
# plausible.
#
# The condition is therefore not "no offset" but "the physical displacements lie on a lattice",
# which is exactly what this checks. It admits the two cases that are actually correct:
#
#   * a single tree, where `offset` is zero and every displacement is a lattice vector;
#   * a block tree whose roots are displaced by a MULTIPLE of the halfsize, where the offset folds
#     into the index and canonicalization then acts on the true physical direction.
#
# and refuses the one that is not. Note this is checked per interaction rather than once on
# `offset`: the halfsize halves each level, so an offset that is a lattice vector at a coarse level
# stays one at every finer level, but not the reverse.
"""
    NonLatticeTranslationError <: Exception

Thrown when `SymmetryDirectionInvariancePerLevel` meets a displacement that is not a lattice
vector, which is the one configuration where a lattice symmetry cannot be applied.

A DISTINCT TYPE rather than a plain `error`, because consumers are expected to CATCH this
specifically and fall back to `DirectionInvariancePerLevel` (`MLFMA` does). A broad `catch` around
a collection construction would swallow method errors, out-of-memory and typos as "not lattice
aligned" and silently produce an unreduced collection with no indication why, which is the failure
mode this check exists to prevent.

PARAMETERIZED ON THE TREE'S OWN COORDINATE TYPE, not fixed to `Float64`. The fields were
`Vector{Float64}`/`Float64`, which did NOT break the catchable-type contract, because Julia's default
constructor converts, so a `Float32`, `BigFloat` or `Rational` tree still produced a genuine,
catchable `NonLatticeTranslationError`. What it broke was the MESSAGE: a `Float32` tree's halfsize
of `0.8f0` was widened and printed as `0.800000011920929`, showing the reader the binary expansion
of a number their input never contained. `BigFloat` lost precision outright, and a coordinate type
not convertible to `Float64` at all would have raised a conversion error in place of this one.
"""
struct NonLatticeTranslationError{V<:AbstractVector{<:Real},T<:Real} <: Exception
    receivingnode::Int
    translatingnode::Int
    translation::V
    halfsize::T
    coordinate::Int
    # Left untyped: the offset is a tree-center difference whose type follows the two trees, and
    # nothing here reads it beyond interpolating it into the message.
    offset::Any
end

function Base.showerror(io::IO, e::NonLatticeTranslationError)
    return print(
        io,
        "SymmetryDirectionInvariancePerLevel: the displacement between nodes " *
        "$(e.receivingnode) and $(e.translatingnode) is $(e.translation), which is not a " *
        "multiple of the level halfsize $(e.halfsize) in coordinate $(e.coordinate) " *
        "($(e.translation[e.coordinate] / e.halfsize) halfsizes). The two trees' roots are " *
        "offset by $(e.offset), and a symmetry cannot be applied across an offset that is not " *
        "itself a lattice vector: the stored direction would be transformed to " *
        "`S*canonical + S*offset` where the interaction needs `S*canonical + offset`. Use " *
        "`DirectionInvariancePerLevel`, which deduplicates by displacement alone and is correct " *
        "for any offset.",
    )
end

function _symmetrylatticeindex!(
    index, translation, halfsize, offset, receivingnode, translatingnode
)
    for i in eachindex(translation)
        scaled = translation[i] / halfsize
        index[i] = round(Int, scaled)
        isapprox(scaled, index[i]; atol=1e-6) || throw(
            NonLatticeTranslationError(
                receivingnode,
                translatingnode,
                collect(translation),
                halfsize,
                i,
                offset,
            ),
        )
    end
    return index
end

"""
    _translations(..., ::SymmetryDirectionInvariancePerLevel; offset=zero(ELTYPE))

Deduplicate `TwoNTree` translations by their lattice displacement up to a symmetry orbit.

The key is the CANONICAL displacement rather than the displacement itself, and the symmetry
mapping the canonical back to the actual is recorded per interaction, so
`q == applysymmetry(group[symmetryID], canonical)` holds for every interaction.

Unlike [`DirectionInvariancePerLevel`](@ref), this requires every displacement to be a lattice
vector. See `_symmetrylatticeindex!` for why a block tree with arbitrarily offset roots is
refused rather than silently reduced.
"""
function _translations(
    translatingplan::AbstractPlan,
    relevantlevels,
    receivecenter,
    translatingcenter,
    receivinghalfsize,
    receivinglevel,
    ::Val{ELTYPE},
    translationtrait::SymmetryDirectionInvariancePerLevel;
    offset=zero(ELTYPE),
) where {ELTYPE}
    D = length(ELTYPE)
    group = symmetrygroup(translationtrait.policy, Val(D))

    translationinfos = [
        @NamedTuple{
            receivingnode::Int, translatingnode::Int, translationID::Int, symmetryID::Int
        }[] for _ in relevantlevels
    ]

    translationID = 1
    temptranslationID = 0

    # keys: [canonicaltranslationperhalfsize..., level]
    translationsdict = Dict{Vector{Int},Int}()

    translation = Vector{eltype(ELTYPE)}(undef, D)
    translationperhalfsize = Vector{Int}(undef, D)
    canonical = Vector{Int}(undef, D)
    canonicalscratch = Vector{Int}(undef, D)
    canonicallevel = Vector{Int}(undef, D + 1)
    vectorindices = 1:D

    halfsizeperlevel = Vector{eltype(ELTYPE)}(undef, length(relevantlevels))

    foreachtranslationpair(
        translatingplan, relevantlevels, receivinglevel
    ) do rlevelid, receivingnode, translatingnode
        rhalfsize = receivinghalfsize(receivingnode)

        halfsizeperlevel[rlevelid] = rhalfsize
        canonicallevel[end] = rlevelid

        # NOT `- offset`. The canonicalization below acts on this index, so it has to be the index
        # of the TRUE physical displacement; an offset folded out here would be reintroduced
        # untransformed. Where the offset is a lattice vector it is already inside this index.
        translation .= receivecenter(receivingnode) - translatingcenter(translatingnode)

        _symmetrylatticeindex!(
            translationperhalfsize,
            translation,
            rhalfsize,
            offset,
            receivingnode,
            translatingnode,
        )

        symmetryID = canonicalizetranslation!(
            canonical, canonicalscratch, translationperhalfsize, group
        )
        view(canonicallevel, vectorindices) .= canonical

        if haskey(translationsdict, canonicallevel)
            temptranslationID = translationsdict[canonicallevel]
        else
            translationsdict[deepcopy(canonicallevel)] = translationID
            temptranslationID = translationID
            translationID += 1
        end

        appendsymmetrytranslationinfo!(
            translationinfos,
            rlevelid,
            receivingnode,
            translatingnode,
            temptranslationID,
            symmetryID,
        )
        return nothing
    end

    ntranslations = length(keys(translationsdict))
    translationarray = Vector{ELTYPE}(undef, ntranslations)

    # NO `+ offset`, matching the index above: the key already is the true physical displacement in
    # halfsize units, so adding the offset back would double-count it.
    for key in keys(translationsdict)
        translationarray[translationsdict[key]] = ELTYPE(
            halfsizeperlevel[key[end]] .* view(key, vectorindices)
        )
    end

    return translationinfos, translationarray
end

function _translations(
    translatingplan::AbstractPlan,
    relevantlevels,
    receivecenter,
    translatingcenter,
    receivinghalfsize,
    receivinglevel,
    ::Val{ELTYPE},
    ::DirectionInvariancePerLevel;
    offset=zero(ELTYPE),
) where {ELTYPE}
    translationinfos = [
        @NamedTuple{receivingnode::Int, translatingnode::Int, translationID::Int}[] for
        _ in relevantlevels
    ]

    translationID = 1
    temptranslationID = 0

    # translations are stored as a multiple of the halfsize at the particular level
    # in the case of two distinct trees an offset has to be subtracted
    # keys: [translationperhalfsize..., level]
    translationsdict = Dict{Vector{Int},Int}()

    # array used as storage for translations
    translation = Vector{eltype(ELTYPE)}(undef, length(ELTYPE))
    # array used as storage for (translation-offset) per halfsize
    translationperhalfsize = Vector{Int}(undef, length(ELTYPE))
    # array used as storage for keys
    translationperhalfsizelevel = Vector{Int}(undef, length(ELTYPE) + 1)
    vectorindices = 1:length(ELTYPE)

    halfsizeperlevel = Vector{eltype(ELTYPE)}(undef, length(relevantlevels))

    foreachtranslationpair(
        translatingplan, relevantlevels, receivinglevel
    ) do rlevelid, receivingnode, translatingnode
        rhalfsize = receivinghalfsize(receivingnode)

        halfsizeperlevel[rlevelid] = rhalfsize
        translationperhalfsizelevel[end] = rlevelid

        translation .=
            receivecenter(receivingnode) - translatingcenter(translatingnode) - offset

        for i in eachindex(translation)
            translationperhalfsize[i] = round(Int, translation[i] / rhalfsize)
        end

        view(translationperhalfsizelevel, vectorindices) .= translationperhalfsize

        if haskey(translationsdict, translationperhalfsizelevel)
            temptranslationID = translationsdict[translationperhalfsizelevel]
        else
            translationsdict[deepcopy(translationperhalfsizelevel)] = translationID
            temptranslationID = translationID
            translationID += 1
        end

        appendtranslationinfo!(
            translationinfos, rlevelid, receivingnode, translatingnode, temptranslationID
        )
        return nothing
    end
    ntranslations = length(keys(translationsdict))
    translationarray = Vector{ELTYPE}(undef, ntranslations)

    for key in keys(translationsdict)
        translationarray[translationsdict[key]] =
            ELTYPE(halfsizeperlevel[key[end]] .* view(key, vectorindices)) + offset
    end

    return translationinfos, translationarray
end

"""
    _translations(..., ::DirectionInvariance; offset=zero(ELTYPE))

Deduplicate `TwoNTree` translations by integer offsets from the minimum
halfsize shared by all relevant levels.
"""
function _translations(
    translatingplan::AbstractPlan,
    relevantlevels,
    receivecenter,
    translatingcenter,
    minhalfsize,
    receivinglevel,
    ::Val{ELTYPE},
    ::DirectionInvariance;
    offset=zero(ELTYPE),
) where {ELTYPE}
    translationinfos = [
        @NamedTuple{receivingnode::Int, translatingnode::Int, translationID::Int}[] for
        _ in relevantlevels
    ]

    translationID = 1
    temptranslationID = 0

    translationsdict = Dict{Vector{Int},Int}()

    # array used as storage for translations
    translation = Vector{eltype(ELTYPE)}(undef, length(ELTYPE))
    # array used as storage for (translation-offset) per halfsize
    translationperhalfsize = Vector{Int}(undef, length(ELTYPE))

    foreachtranslationpair(
        translatingplan, relevantlevels, receivinglevel
    ) do rlevelid, receivingnode, translatingnode
        translation .=
            receivecenter(receivingnode) - translatingcenter(translatingnode) - offset

        for i in eachindex(translation)
            translationperhalfsize[i] = round(Int, translation[i] / minhalfsize)
        end

        if haskey(translationsdict, translationperhalfsize)
            temptranslationID = translationsdict[translationperhalfsize]
        else
            translationsdict[deepcopy(translationperhalfsize)] = translationID
            temptranslationID = translationID
            translationID += 1
        end

        appendtranslationinfo!(
            translationinfos, rlevelid, receivingnode, translatingnode, temptranslationID
        )
        return nothing
    end

    ntranslations = length(keys(translationsdict))
    translationarray = Vector{ELTYPE}(undef, ntranslations)

    for key in keys(translationsdict)
        translationarray[translationsdict[key]] = ELTYPE(minhalfsize .* key) + offset
    end

    return translationinfos, translationarray
end
