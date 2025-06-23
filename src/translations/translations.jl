struct CenterFunctor{T}
    tree::T
end

function (f::CenterFunctor)(node::Int)
    return H2Trees.center(f.tree, node)
end

struct HalfSizeFunctor{T}
    tree::T
end

function (f::HalfSizeFunctor)(node::Int)
    return H2Trees.halfsize(f.tree, node)
end

struct LevelFunctor{T}
    tree::T
end

function (f::LevelFunctor)(node::Int)
    return H2Trees.level(f.tree, node)
end

function translations(tree, translatingplan::AbstractPlan, translationtrait)
    @assert istranslatingplan(translatingplan)
    return translations(tree, treetrait(tree), translatingplan, translationtrait)
end

function translations(
    tree, ::AbstractTreeTrait, translatingplan::AbstractPlan, translationtrait
)
    relevantlevels = mintranslationlevel(translatingplan):(H2Trees.levels(tree)[end])

    return _translations(
        translatingplan,
        relevantlevels,
        CenterFunctor(tree),
        CenterFunctor(tree),
        LevelFunctor(tree),
        Val{eltype(tree)}(),
        translationtrait,
    )
end

function translations(
    tree, ::isTwoNTree, translatingplan::AbstractPlan, ::DirectionInvariancePerLevel
)
    relevantlevels = mintranslationlevel(translatingplan):(levels(tree)[end])

    return _translations(
        translatingplan,
        relevantlevels,
        CenterFunctor(tree),
        CenterFunctor(tree),
        HalfSizeFunctor(tree),
        LevelFunctor(tree),
        Val{eltype(tree)}(),
        DirectionInvariancePerLevel(),
    )
end

function translations(
    tree, ::isTwoNTree, translatingplan::AbstractPlan, ::DirectionInvariance
)
    relevantlevels = mintranslationlevel(translatingplan):(levels(tree)[end])

    return _translations(
        translatingplan,
        relevantlevels,
        CenterFunctor(tree),
        CenterFunctor(tree),
        H2Trees.minhalfsize(tree),
        LevelFunctor(tree),
        Val{eltype(tree)}(),
        DirectionInvariance(),
    )
end

function translations(tree, ::isBlockTree, translatingplan::AbstractPlan, translationtrait)
    return translations(
        testtree(tree),
        trialtree(tree),
        treetrait(testtree(tree)),
        treetrait(trialtree(tree)),
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

    return _translations(
        translatingplan,
        relevantlevels,
        CenterFunctor(receivetree),
        CenterFunctor(translatingtree),
        LevelFunctor(receivetree),
        Val{promote_type(eltype(receivetree), eltype(translatingtree))}(),
        translationtrait,
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

    return _translations(
        translatingplan,
        relevantlevels,
        CenterFunctor(receivetree),
        CenterFunctor(translatingtree),
        HalfSizeFunctor(receivetree),
        LevelFunctor(receivetree),
        Val{promote_type(eltype(receivetree), eltype(translatingtree))}(),
        DirectionInvariancePerLevel();
        offset=H2Trees.center(receivetree, H2Trees.root(receivetree)) -
               H2Trees.center(translatingtree, H2Trees.root(translatingtree)),
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

    return _translations(
        translatingplan,
        relevantlevels,
        CenterFunctor(receivetree),
        CenterFunctor(translatingtree),
        min(H2Trees.minhalfsize(receivetree), H2Trees.minhalfsize(translatingtree)),
        LevelFunctor(receivetree),
        Val{promote_type(eltype(receivetree), eltype(translatingtree))}(),
        DirectionInvariance();
        offset=H2Trees.center(receivetree, H2Trees.root(receivetree)) -
               H2Trees.center(translatingtree, H2Trees.root(translatingtree)),
    )
end

function _translations(
    translatingplan::AbstractPlan,
    relevantlevels,
    receivecenter,
    translatingcenter,
    receivinglevel,
    ::Val{ELTYPE},
    ::AllTranslations,
) where {ELTYPE}
    relevantlevelsdict = Dict(zip(relevantlevels, collect(eachindex(relevantlevels))))

    # count number of translations on each level
    ntranslationsperlevel = zeros(Int, length(relevantlevels))
    for level in levels(translatingplan)
        for receivingnode in receivingnodes(translatingplan, level)
            ntranslationsperlevel[relevantlevelsdict[receivinglevel(receivingnode)]] += length(
                translatingplan[receivingnode, level]
            )
        end
    end

    # allocate required memory
    translationinfos = Vector{
        Vector{@NamedTuple{receivingnode::Int,translatingnode::Int,translationID::Int}}
    }(
        undef, length(relevantlevels)
    )
    for level in relevantlevels
        translationinfos[relevantlevelsdict[level]] = Vector{
            @NamedTuple{receivingnode::Int,translatingnode::Int,translationID::Int}
        }(
            undef, ntranslationsperlevel[relevantlevelsdict[level]]
        )
    end
    translations = Vector{ELTYPE}(undef, sum(ntranslationsperlevel))

    # compute translations
    translationID = 1
    translationIDlevel = ones(Int, length(relevantlevels))
    for level in relevantlevels
        for receivingnode in receivingnodes(translatingplan, level)
            rlevelid = relevantlevelsdict[receivinglevel(receivingnode)]
            rcenter = receivecenter(receivingnode)
            for translationnode in translatingplan[receivingnode, level]
                translation = rcenter - translatingcenter(translationnode)

                translations[translationID] = translation

                translationinfo = (
                    receivingnode=receivingnode,
                    translatingnode=translationnode,
                    translationID=translationID,
                )
                translationinfos[rlevelid][translationIDlevel[rlevelid]] = translationinfo

                translationID += 1
                translationIDlevel[rlevelid] += 1
            end
        end
    end

    return translationinfos, translations
end

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
    relevantlevelsdict = Dict(zip(relevantlevels, collect(eachindex(relevantlevels))))

    translations = [ELTYPE[] for _ in relevantlevels]
    translationIDs = [Int[] for _ in relevantlevels]
    translationinfos = [
        @NamedTuple{receivingnode::Int, translatingnode::Int, translationID::Int}[] for
        _ in relevantlevels
    ]

    translationID = 1
    temptranslationID = 0

    for level in relevantlevels
        for receivenode in receivingnodes(translatingplan, level)
            rlevelid = relevantlevelsdict[receivinglevel(receivenode)]
            rcenter = receivecenter(receivenode)
            for translationnode in translatingplan[receivenode, level]
                translation = rcenter - translatingcenter(translationnode)
                translationindexinarray = findfirst(
                    isapprox(translation), translations[rlevelid]
                )

                if isnothing(translationindexinarray)
                    push!(translations[rlevelid], translation)
                    push!(translationIDs[rlevelid], translationID)
                    temptranslationID = translationID
                    translationID += 1
                else
                    temptranslationID = translationIDs[rlevelid][translationindexinarray]
                end
                translationinfo = (
                    receivingnode=receivenode,
                    translatingnode=translationnode,
                    translationID=temptranslationID,
                )

                push!(translationinfos[rlevelid], translationinfo)
            end
        end
    end

    translationarray = Vector{ELTYPE}(undef, translationID - 1)

    for i in eachindex(translations)
        for j in eachindex(translations[i])
            translationarray[translationIDs[i][j]] = translations[i][j]
        end
    end

    return translationinfos, translationarray
end

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
    relevantlevelsdict = Dict(zip(relevantlevels, collect(eachindex(relevantlevels))))

    translations = ELTYPE[]
    translationIDs = Int[]
    translationinfos = [
        @NamedTuple{receivingnode::Int, translatingnode::Int, translationID::Int}[] for
        _ in relevantlevels
    ]

    translationID = 1
    temptranslationID = 0

    for level in relevantlevels
        for receivenode in receivingnodes(translatingplan, level)
            rlevelid = relevantlevelsdict[receivinglevel(receivenode)]
            rcenter = receivecenter(receivenode)
            for translationnode in translatingplan[receivenode, level]
                translation = rcenter - translatingcenter(translationnode)
                translationindexinarray = findfirst(isapprox(translation), translations)

                if isnothing(translationindexinarray)
                    push!(translations, translation)
                    push!(translationIDs, translationID)
                    temptranslationID = translationID
                    translationID += 1
                else
                    temptranslationID = translationIDs[translationindexinarray]
                end

                translationinfo = (
                    receivingnode=receivenode,
                    translatingnode=translationnode,
                    translationID=temptranslationID,
                )

                push!(translationinfos[rlevelid], translationinfo)
            end
        end
    end
    return translationinfos, translations
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
    relevantlevelsdict = Dict(zip(relevantlevels, collect(eachindex(relevantlevels))))

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

    for level in relevantlevels
        for receivenode in receivingnodes(translatingplan, level)
            rlevelid = relevantlevelsdict[receivinglevel(receivenode)]
            rcenter = receivecenter(receivenode)
            rhalfsize = receivinghalfsize(receivenode)

            halfsizeperlevel[rlevelid] = rhalfsize
            translationperhalfsizelevel[end] = rlevelid

            for translationnode in translatingplan[receivenode, level]
                translation .= rcenter - translatingcenter(translationnode) - offset

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

                translationinfo = (
                    receivingnode=receivenode,
                    translatingnode=translationnode,
                    translationID=temptranslationID,
                )

                push!(translationinfos[rlevelid], translationinfo)
            end
        end
    end
    ntranslations = length(keys(translationsdict))
    translationarray = Vector{ELTYPE}(undef, ntranslations)

    for key in keys(translationsdict)
        translationarray[translationsdict[key]] =
            ELTYPE(halfsizeperlevel[key[end]] .* view(key, vectorindices)) + offset
    end

    return translationinfos, translationarray
end

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
    relevantlevelsdict = Dict(zip(relevantlevels, collect(eachindex(relevantlevels))))

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

    for level in relevantlevels
        for receivenode in receivingnodes(translatingplan, level)
            rlevelid = relevantlevelsdict[receivinglevel(receivenode)]
            rcenter = receivecenter(receivenode)

            for translationnode in translatingplan[receivenode, level]
                translation .= rcenter - translatingcenter(translationnode) - offset

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

                translationinfo = (
                    receivingnode=receivenode,
                    translatingnode=translationnode,
                    translationID=temptranslationID,
                )

                push!(translationinfos[rlevelid], translationinfo)
            end
        end
    end

    ntranslations = length(keys(translationsdict))
    translationarray = Vector{ELTYPE}(undef, ntranslations)

    for key in keys(translationsdict)
        translationarray[translationsdict[key]] = ELTYPE(minhalfsize .* key) + offset
    end

    return translationinfos, translationarray
end
