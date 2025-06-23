module H2Trees
using StaticArrays
using LinearAlgebra
import Base.Threads: @threads
using ChunkSplitters
using ProgressMeter
using Match

function QuadPointsTree end # requires BEAST to load
export QuadPointsTree

abstract type H2ClusterTree end

function traceball end # requires PlotlyJS to load
function tracecube end # requires PlotlyJS to load

export traceball, tracecube

function progressbar(workload, verbose; kwargs...)
    return Progress(
        workload; barglyphs=BarGlyphs("[=> ]"), color=:white, enabled=verbose, kwargs...
    )
end

include("treetraits.jl")

abstract type FarMulMode end

"""
    AggregateMode <: FarMulMode

This mode uses `AggregatePlan` and `DisaggregateTranslatePlan` to perform the
farmultiplication.
"""
struct AggregateMode <: FarMulMode end

"""
    AggregateTranslateMode <: FarMulMode

This mode uses `AggregateTranslatePlan` and `DisaggregatePlan` to perform the
farmultiplication.
"""
struct AggregateTranslateMode <: FarMulMode end

function Base.adjoint(::AggregateMode)
    return AggregateTranslateMode()
end

function Base.adjoint(::AggregateTranslateMode)
    return AggregateMode()
end

export treetrait
export isTwoNTree
export isBlockTree
export UniquePoints
export NonUniquePoints

struct Node{D}
    data::D
    next_sibling::Int
    parent::Int
    first_child::Int
end

struct BoxData{N,T}
    sector::Int
    values::Vector{Int}
    center::SVector{N,T}
    halfsize::T
    level::Int
end

struct ParametricBoundingBallData{N,T}
    values::Vector{Int}
    center::SVector{N,T}
    radius::T
    level::Int
    parametricsector::Int
    parametricnode::Int
    patchID::Int
end

include("bounding/boundingbox.jl")

include("iterators/NodeFilterIterator.jl")

include("iterators/nearfar/NearNodeIterator.jl")
include("iterators/nearfar/nearinteractions.jl")
# include("iterators/nearfar/corrector.jl")
include("iterators/nearfar/isnear.jl")

include("iterators/WellSeparatedIterator.jl")
include("iterators/AllLeavesTranslationsIterator.jl")

include("plans/plans.jl")

include("translations/translationtraits.jl")
include("translations/translations.jl")

include("printing.jl")

include("testingutils/testingutils.jl")

export AggregatePlan, AggregateTranslatePlan, DisaggregatePlan, DisaggregateTranslatePlan

function leafclusters(tree)
    clusters = Vector{Vector{Int}}(undef, length(H2Trees.leaves(tree)))
    for (i, leaf) in enumerate(H2Trees.leaves(tree))
        clusters[i] = H2Trees.values(tree, leaf)
    end

    return clusters
end

# API ######################################################################################

export ParametricBoundingBallData
export BoxData
export Node

function data(tree, node::Int)
    return tree(node).data
end

"""
    values(tree, node::Int)

Returns the values stored in the given `node` of the `tree`. If the `node` is a leaf node,
it returns the values directly. Otherwise, it recursively collects the values from all the
leaf nodes in the subtree rooted at the given `node`.

# Arguments

  - `tree`: The H2 tree.
  - `node::Int`: The index of the node.

# Returns

An array of values stored in the given `node` or its subtree.
"""
function values(tree, node::Int)
    iszero(firstchild(tree, node)) && return H2Trees.values(data(tree, node))

    values = Int[]
    for i in H2Trees.leaves(tree, node)
        append!(values, H2Trees.values(tree, Int(i)))
    end
    return values
end

function values(data::Union{BoxData,ParametricBoundingBallData})
    return data.values
end

function root(tree)
    return tree.root
end

function level(tree, nodeid::Int)
    return level(tree(nodeid).data) #TODO: look at code duplications: level, sector, data: solve with metaprogramming
end

function level(data::Union{BoxData,ParametricBoundingBallData})
    return data.level
end

"""
    minimumlevel(tree)

Get the minimum level of a tree, which is the level of the root node. This is not
necessarily the level 1.

# Arguments

  - `tree`: The tree.

# Returns

The minimum level of the tree.
"""
function minimumlevel(tree)
    return level(tree(root(tree)).data)
end

function levels(tree)
    return (1:length(nodesatlevel(tree))) .+ (H2Trees.minimumlevel(tree) - 1)
end

function numberoflevels(tree)
    return length(levels(tree))
end

function parent(tree, node::Int)
    return tree(node).parent
end

function nextsibling(tree, node::Int)
    return tree(node).next_sibling
end

function firstchild(tree, node::Int)
    return tree(node).first_child
end

function numberofnodes(tree)
    return length(tree.nodes)
end
function leaves(tree, node::Int) end

function isleaf(tree, node::Int)
    return iszero(tree(node).first_child)
end
# returns the leaf node that contains the given point
"""
    findleafnode(tree, value::Int)

Find the leaf node in the given `tree` that contains the specified `value`.

# Arguments

  - `tree`: The tree to search in.
  - `value`: The value to search for.

# Returns

  - The leaf node that contains the `value`, or `0` if not found.
"""
function findleafnode(tree, value::Int)
    for leaf in H2Trees.leaves(tree)
        (value ∈ H2Trees.values(tree, leaf)) && return leaf
    end
    return 0
end

function children(tree, node::Int=root(tree))
    return ChildIterator(tree, node)
end

function center(tree, nodeid::Int=root(tree))
    return center(tree(nodeid).data)
end

function center(data::Union{BoxData,ParametricBoundingBallData})
    return data.center
end

function sector(tree, node::Int)
    return sector(tree(node).data)
end

function sector(data::BoxData)
    return data.sector
end

function sector(data::ParametricBoundingBallData)
    return data.parametricsector
end

function halfsize(tree, nodeid::Int=root(tree))
    return halfsize(data(tree, nodeid))
end

function halfsize(data::BoxData)
    return data.halfsize
end

function halfsizes(tree)
    halfsizes = eltype(eltype(tree))[]
    for level in H2Trees.levels(tree)
        for node in H2Trees.nodesatlevel(tree, level)
            push!(halfsizes, H2Trees.halfsize(tree, node))
            break
        end
    end
    return halfsizes
end

function minhalfsize(tree)
    nodesatlowestlevel = nodesatlevel(tree)[end]
    nodesatlowestlevel == Int[] && return 0

    return data(tree, nodesatlowestlevel[begin]).halfsize
end

function radius(tree, node::Int)
    return radius(tree(node).data)
end

function radius(data::ParametricBoundingBallData)
    return data.radius
end

function treewithmorelevels(tree)
    if length(levels(trialtree(tree))) >= length(levels(testtree(tree)))
        return trialtree(tree)
    else
        return testtree(tree)
    end
end

function parametrictree(tree)
    return tree.parametrictree
end

function parametricnode(tree, node)
    return parametricnode(tree(node).data)
end

function parametricnode(data::ParametricBoundingBallData)
    return data.parametricnode
end

function patchID(tree, node)
    return patchID(tree(node).data)
end

function patchID(data::ParametricBoundingBallData)
    return data.patchID
end

function samelevelnodes(tree, nodeid::Int)
    level = H2Trees.level(tree, nodeid)
    return nodesatlevel(tree, level)
end

function nodesatlevel(tree)
    return tree.nodesatlevel
end

function nodesatlevel(tree, level::Int)
    levels = H2Trees.levels(tree)
    level < levels[begin] && return Int[]
    level > levels[end] && return Int[]
    return H2Trees.nodesatlevel(tree)[leveltolevelid(tree, level)]
end

"""
    LevelIterator(tree, level::Int)

Return an iterator over the nodes at the specified `level` in the `tree`.

# Arguments

  - `tree`: The tree object.
  - `level`: The level at which to iterate.

# Returns

An iterator over the nodes at the specified level.
"""
function LevelIterator(tree, level::Int)
    return nodesatlevel(tree, level)
end

"""
    SameLevelIterator(tree, node::Int)

Returns an iterator over the nodes at the same level as `node` in the `tree`.

# Arguments

  - `tree`: The tree structure.
  - `node`: The node for which to find the same level nodes.

# Returns

An iterator over the nodes at the same level as `node`.
"""
function SameLevelIterator(tree, node::Int)
    return samelevelnodes(tree, node)
end

function _adjustnodesatlevels!(tree)
    empty!(tree.nodesatlevel)
    for node in H2Trees.DepthFirstIterator(tree, H2Trees.root(tree))
        nodelevel = H2Trees.level(tree, node)
        minlevel = H2Trees.minimumlevel(tree)
        if nodelevel - minlevel + 1 > length(nodesatlevel(tree))
            numberofmissinglevels = nodelevel - minlevel - length(nodesatlevel(tree))
            for _ in 1:numberofmissinglevels
                push!(tree.nodesatlevel, Int[])
            end
            push!(tree.nodesatlevel, [node])

        else
            append!(tree.nodesatlevel[leveltolevelid(tree, nodelevel)], node)
        end
    end

    for nodesatlevel in tree.nodesatlevel
        sort!(nodesatlevel)
    end
end

function numberofvalues(tree)
    maxvalue = 0
    for leaf in H2Trees.leaves(tree)
        maxvalue = max(maxvalue, maximum(H2Trees.values(tree, leaf)))
    end

    return maxvalue
end

function valuesatnodes(tree; numberofvalues=H2Trees.numberofvalues(tree))
    nodes = Vector{Vector{Int}}(undef, numberofvalues)

    for leaf in H2Trees.leaves(tree)
        for value in H2Trees.values(tree, leaf)
            if isassigned(nodes, value)
                push!(nodes[value], leaf)
            else
                nodes[value] = [leaf]
            end
        end
    end

    for value in 1:numberofvalues
        if isassigned(nodes, value)
            sort!(nodes[value])
        else
            nodes[value] = Int[]
        end
    end
    return nodes
end

function nodesatvalues(tree, boxes=H2Trees.valuesatnodes(tree))
    boxesdict = Dict{Vector{Int},Vector{Int}}()
    for (value, box) in enumerate(boxes)
        if haskey(boxesdict, box)
            push!(boxesdict[box], value)
        else
            boxesdict[box] = [value]
        end
    end

    return boxesdict
end

function uniquepointstree(tree)
    return uniquepointstree(tree, treetrait(tree))
end

function uniquepointstree(tree, ::isBlockTree)
    return uniquepointstree(testtree(tree)) && uniquepointstree(trialtree(tree))
end

function uniquepointstree(tree, ::AbstractTreeTrait)
    return true
end

"""
    leveltolevelid(tree, level::Int)

Converts a level in the tree to its corresponding level ID. This is relevant since the first
level might not be level 1.

# Arguments

  - `tree`: The tree object.
  - `level`: The level to convert.

# Returns

The level ID corresponding to the given level.
"""
function leveltolevelid(tree, level::Int)
    return level - H2Trees.minimumlevel(tree) + 1
end

function levelindex(tree, node::Int)
    return leveltolevelid(tree, H2Trees.level(tree, node))
end

#TODO: rethink this naming
function checkbalancedtree(tree)
    leaflevel = level(tree, H2Trees.leaves(tree)[1])
    for node in H2Trees.leaves(tree)
        leaflevel != level(tree, node) && return false
    end
    return true
end

function computevectorbuffers(tree, T)
    return computevectorbuffers(tree, treetrait(tree), T)
end

function computevectorbuffers(tree, ::isBlockTree, T)
    return computevectorbuffers(testtree(tree), T), computevectorbuffers(trialtree(tree), T)
end

function computevectorbuffers(tree, ::Any, ::Type{T}) where {T}
    vectors = Vector{Vector{T}}(undef, length(H2Trees.leaves(tree)))

    for (i, leaf) in enumerate(H2Trees.leaves(tree))
        vectors[i] = Vector{T}(undef, length(H2Trees.values(tree, leaf)))
    end
    return Dict(zip(H2Trees.leaves(tree), vectors))
end

function isuppertreelevel(tree, level::Int)
    return level ≤ hybridlevel(tree)
end

function islowertreelevel(tree, level::Int)
    return !isuppertreelevel(tree, level)
end

function isuppertreenode(tree, node::Int)
    return isuppertreelevel(tree, level(tree, node))
end

function islowertreenode(tree, node::Int)
    return !H2Trees.isuppertreenode(tree, node)
end

include("trees/clustertrees.jl")
include("trees/TwoNTree.jl")
include("trees/SimpleHybridTree.jl")
include("trees/BoundingBallTree.jl")
include("trees/BlockTree.jl")

export TwoNTree, BlockTree
end
