struct NoCorrectionPlan end

struct CorrectionPlan{D} <: AbstractPlan
    receivingnodes::Dict{Int,Vector{Int}}
    correctiondict::D
    leaflevel::Int

    function CorrectionPlan(tree; kwargs...)
        H2Trees.arepointsunique(tree) && return NoCorrectionPlan()
        cinfo = H2Trees.correctioninformation(tree; kwargs...)

        receivingnodes = Dict{Int,Vector{Int}}()
        for (translatingnode, receivingnode) in keys(cinfo)
            !haskey(receivingnodes, receivingnode) &&
                (receivingnodes[receivingnode] = Int[])
            push!(receivingnodes[receivingnode], translatingnode)
        end

        return new{typeof(cinfo)}(receivingnodes, cinfo, H2Trees.levels(tree)[end])
    end
end

function levels(plan::CorrectionPlan)
    return (plan.leaflevel):(plan.leaflevel)
end

function plantranslationtrait(::CorrectionPlan)
    return IsTranslatingPlan()
end

function mintranslationlevel(plan::CorrectionPlan)
    return plan.leaflevel
end

function receivingnodes(plan::CorrectionPlan, level)
    !(level == plan.leaflevel) && return Int[]

    return collect(keys(plan.receivingnodes))
end

function Base.getindex(plan::CorrectionPlan, receivingnode::Int, level::Int)
    !(level == plan.leaflevel) && return Int[]

    tfnodes = plan.receivingnodes

    if haskey(tfnodes, receivingnode)
        return tfnodes[receivingnode]
    else
        return Int[]
    end
end
