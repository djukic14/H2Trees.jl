
struct CorrectNearInteractionsFunctor{N,S}
    nearinteractions::N
    symmetric::S
end

function (f::CorrectNearInteractionsFunctor)(testvalue, trialvalues, correction)
    return correct(
        f.nearinteractions, testvalue, trialvalues, correction; symmetric=f.symmetric
    )
end

function correct(nearinteractions, testvalue, trialvalues, correction; symmetric=false)
    if symmetric
        nearinteractions[trialvalues, testvalue] -= correction
        nearinteractions[testvalue, trialvalues] -= correction
    else
        nearinteractions[trialvalues, testvalue] -= correction
    end
end

struct NoNearCorrection end

struct NearCorrection{T}
    tree::T
end

function correct(tree, T)
    return !(corrector(tree, T) == NoNearCorrection())
end

function corrector(tree, ::Any)
    return corrector(tree)
end

function corrector(tree)
    return treecorrector(tree, uniquepointstreetrait(tree))
end

function treecorrector(::Any, ::UniquePoints)
    return NoNearCorrection()
end

function treecorrector(tree, ::NonUniquePoints)
    return NearCorrection(tree)
end

function (::NoNearCorrection)(::Any, ::Any; kwargs...) end

function (f::NearCorrection)(blockassembler, corrector; isnear=isnear, verbose=true)
    return nearinteractionscorrection!(
        f.tree, blockassembler, corrector; isnear=isnear, verbose=verbose
    )
end

function correctioninformation(
    tree;
    isnear=isnear,
    valuesatnodes=H2Trees.valuesatnodes(tree),
    nodesatvalues=H2Trees.nodesatvalues(tree, valuesatnodes),
)
    lk = Threads.SpinLock()

    correctioninformation = Dict{NTuple{2,Int},Dict{Vector{Int},Vector{Int}}}()

    @threads :static for (boxes, sharedvalues) in collect(nodesatvalues)
        nearvalues = Set{Int}()
        for box in boxes
            for nearnode in NearNodeIterator(tree, box; isnear=isnear)
                setappend!(nearvalues, values(tree, nearnode))
            end
        end

        extendednearboxes = Set{Int}()
        for nearvalue in nearvalues
            setappend!(extendednearboxes, valuesatnodes[nearvalue])
        end

        for box in boxes
            for otherbox in extendednearboxes
                isnear(tree, box, otherbox) && continue

                otherboxvalues = values(tree, otherbox)
                nearotherboxvalues = Int[]
                for val in otherboxvalues
                    val ∉ nearvalues && continue
                    push!(nearotherboxvalues, val)
                end

                isempty(nearotherboxvalues) && continue

                lock(lk) do
                    if haskey(correctioninformation, (box, otherbox))
                        correctioninformation[(box, otherbox)][sharedvalues] =
                            nearotherboxvalues
                    else
                        correctioninformation[(box, otherbox)] = Dict(
                            zip([sharedvalues], [nearotherboxvalues])
                        )
                    end
                end
            end
        end
    end
    return correctioninformation
end

function nearinteractionscorrection!(
    tree,
    blockassembler,
    corrector;
    isnear=isnear,
    verbose=true,
    correctionplan=CorrectionPlan(tree; isnear=isnear),
    symmetric=corrector.symmetric,
)
    lk = Threads.SpinLock()

    cinformation = correctionplan.correctiondict

    combinations = collect(keys(cinformation))

    @threads :static for key in combinations
        correctioninfo = cinformation[key]
        box, otherbox = key
        (symmetric && (box > otherbox)) && continue
        # compute block matrix between box and otherbox

        testvalues = Set{Int}()
        trialvalues = Set{Int}()
        for (sharedvalues, nearotherboxvalues) in correctioninfo
            setappend!(testvalues, sharedvalues)
            setappend!(trialvalues, nearotherboxvalues)
        end

        testvalues = collect(testvalues)
        trialvalues = collect(trialvalues)

        localtestvalues = Int[]
        localtrialvalues = Int[]
        for testvalue in testvalues
            push!(localtestvalues, findfirst(isequal(testvalue), H2Trees.values(tree, box)))
        end

        for trialvalue in trialvalues
            push!(
                localtrialvalues,
                findfirst(isequal(trialvalue), H2Trees.values(tree, otherbox)),
            )
        end

        testdict = Dict(zip(testvalues, 1:length(localtestvalues)))
        trialdict = Dict(zip(trialvalues, 1:length(localtrialvalues)))

        correction = blockassembler(box, otherbox, localtestvalues, localtrialvalues)

        for (sharedvalues, nearotherboxvalues) in correctioninfo
            for val in sharedvalues
                localval = testdict[val]
                for otherval in nearotherboxvalues
                    localotherval = trialdict[otherval]
                    #TODO: this can probably be sorted such that we do not need a lock
                    # for Blocksparse matrices
                    lock(lk) do
                        return corrector(val, otherval, correction[localval, localotherval])
                    end
                end
            end
        end
    end

    return nothing
end
