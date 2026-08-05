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

"""
    _translations(..., ::DirectionInvariancePerLevel; offset=zero(ELTYPE))

Deduplicate `TwoNTree` translations by integer offsets from the level halfsize.

For a block tree with offset roots, `offset` shifts both sides to the same grid
before the integer direction key is computed.
"""
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
