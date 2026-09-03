module H2Trees
using LinearAlgebra
using Graphs
using BoundingSphere
using StaticArrays
using Random
import Base.Threads: @threads

# Self-contained smallest-enclosing-ball-of-balls solver. Included before any tree source so
# the `BoundingBallTree` adapter can call into it. `SEBB` itself must never depend on tree
# types (dependency direction: SEBB <- H2Trees adapter).
include("SEBB/SEBB.jl")
export SEBB

"""
    H2ClusterTree

Abstract supertype for concrete cluster-tree implementations.
"""
abstract type H2ClusterTree end

"""
    traceball(...)
    tracecube(...)

Extension hooks implemented by `H2PlotlyJSTrees`.
"""
function traceball end # requires PlotlyJS to load
function tracecube end # requires PlotlyJS to load

export traceball, tracecube

include("forests/forest.jl")

include("treetraits.jl")

"""
    FarMulMode

Abstract supertype for far-field multiplication traversal modes.
"""
abstract type FarMulMode end

"""
    AggregateMode <: FarMulMode

Far-field traversal mode using [`AggregatePlan`](@ref) with
[`DisaggregateTranslatePlan`](@ref).
"""
struct AggregateMode <: FarMulMode end

"""
    AggregateTranslateMode <: FarMulMode

Far-field traversal mode using [`AggregateTranslatePlan`](@ref) with
[`DisaggregatePlan`](@ref).
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
export SameLevelFilteredIterator, NearNodesAtAnchorLevel, FarNodesAtAnchorLevel
export appendvalues!, foreachvalue, anyvalue, appendnearnodevalues!, appendfarnodevalues!

"""
    Node(data, nextsibling, parent, firstchild)

Topology record stored by cluster trees.
"""
struct Node{D}
    data::D
    nextsibling::Int
    parent::Int
    firstchild::Int
end

"""
    BoxData(sector, values, center, halfsize, level)

Node payload for [`TwoNTree`](@ref)-style box trees.
"""
struct BoxData{N,T}
    sector::Int
    values::Vector{Int}
    center::SVector{N,T}
    halfsize::T
    level::Int
end

"""
    BoundingBallData(values, center, radius, level)

Node payload for [`BoundingBallTree`](@ref)-style ball trees.
"""
struct BoundingBallData{N,T}
    values::Vector{Int}
    center::SVector{N,T}
    radius::T
    level::Int
end

include("hilbert/HilbertOrdering.jl")

include("bounding/boundingbox.jl")

include("iterators/NodeFilterIterator.jl")

include("iterators/nearfar/NearNodeIterator.jl")
include("iterators/nearfar/nearinteractions.jl")
include("iterators/nearfar/isnear.jl")
include("iterators/nearfar/nearlists.jl")

include("iterators/WellSeparatedIterator.jl")

include("plans/plans.jl")

include("translations/latticesymmetry.jl")
include("translations/translationtraits.jl")
include("translations/translations.jl")

export LatticeSymmetry, SymmetryGroup, symmetrygroup
export NoSymmetry, OppositeSymmetry, AxisPreservingSymmetry, FullLatticeSymmetry
export canonicalizetranslation, canonicalizetranslation!, symmetryorbit
export SymmetryDirectionInvariancePerLevel

include("printing.jl")

include("protrusion.jl")
include("treebuilders.jl")

export AggregatePlan, AggregateTranslatePlan, DisaggregatePlan, DisaggregateTranslatePlan
export checkadmissibility, AdmissibilityReport, AdmissibilityFinding
export PlanBuilder, PlanSet, buildplans
export GalerkinPlanFamily, PetrovPlanFamily, planfamily
export trialaggregationplan, testdisaggregationplan, testaggregationplan
export trialdisaggregationplan, relevantlevels, mintranslationlevel
export ownedtree, receivingtree, translatingtree
export refreshreceivingnodes!
export isgalerkin, ispetrov
export AggregateMode, AggregateTranslateMode

"""
    leafclusters(tree)

Return the leaf value vectors of `tree`.
"""
function leafclusters(tree)
    clusters = Vector{Vector{Int}}(undef, length(leaves(tree)))
    for (i, leaf) in enumerate(leaves(tree))
        clusters[i] = values(data(tree, leaf))
    end

    return clusters
end

# API ######################################################################################

export BoundingBallData
export BoxData
export Node
export balanceleaves!
export AutoMinLevel
export TwoNTreeBuilder, BlockTreeBuilder
export BoundingBallTreeBuilder, KMeansTreeBuilder, MetisTreeBuilder, MetisForestBuilder
export SimpleHybridTreeBuilder
export buildtree, buildforest
export NoProtrusionCheck, AutoProtrusionCheck, ProtrusionCheck
export defaultprotrusioncheck, defaultprotrusionfunctor
export TreeIndex, treeindex, depthfirstnodes, rebuildtreeindex!
export levelprotrusions, protrusionreport

"""
    TreeIndex

Cached topology views for a tree.

`TreeIndex` stores nodes grouped by level, depth-first order, leaves, and the
minimum/maximum levels. `rebuildtreeindex!` replaces the whole cached value in
the tree's `Ref`, so callers never observe a half-updated index. The vectors
inside the index are ordinary mutable arrays; treat them as read-only.
"""
struct TreeIndex
    nodes_by_level::Vector{Vector{Int}}
    depthfirstnodes::Vector{Int}
    leaves::Vector{Int}
    minlevel::Int
    maxlevel::Int
end

"""
    TreeIndex(tree)

Compute cached topology views for `tree`.
"""
function TreeIndex(tree)
    depthfirst = collect(Int, DepthFirstIterator(tree, root(tree)))
    minlevel = isempty(depthfirst) ? 0 : minimum(level.(Ref(tree), depthfirst))
    maxlevel = isempty(depthfirst) ? 0 : maximum(level.(Ref(tree), depthfirst))
    nodes_by_level = [Int[] for _ in minlevel:maxlevel]
    leafnodes = Int[]

    for node in depthfirst
        push!(nodes_by_level[level(tree, node) - minlevel + 1], node)
        isleaf(tree, node) && push!(leafnodes, node)
    end

    for nodes in nodes_by_level
        sort!(nodes)
    end

    return TreeIndex(nodes_by_level, depthfirst, leafnodes, minlevel, maxlevel)
end

"""
    treeindex(tree)

Return the cached [`TreeIndex`](@ref) for `tree`.
"""
treeindex(tree) = tree.index[]

"""
    data(tree, node)

Return the data payload stored on `node`.
"""
function data(tree, node::Int)
    return tree(node).data
end

"""
    foreachvalue(f, tree, node::Int)

Call `f(value)` for every value contained in `node`'s subtree.

For a leaf this visits the values physically stored on that node. For an internal node this walks
the child subtrees directly, avoiding the temporary leaf or value vectors materialized by
`leaves(tree, node)` and `values(tree, node)`.
"""
function foreachvalue(f, tree, node::Int)
    if iszero(firstchild(tree, node))
        for value in values(data(tree, node))
            f(value)
        end
    else
        for child in children(tree, node)
            foreachvalue(f, tree, child)
        end
    end
    return nothing
end

"""
    anyvalue(f, tree, node::Int)

Return `true` when `f(value)` is true for any value contained in `node`'s subtree.
"""
function anyvalue(f, tree, node::Int)
    if iszero(firstchild(tree, node))
        for value in values(data(tree, node))
            f(value) && return true
        end
    else
        for child in children(tree, node)
            anyvalue(f, tree, child) && return true
        end
    end
    return false
end

"""
    appendvalues!(out, tree, node::Int)

Append all values contained in `node`'s subtree to `out`.

For a leaf this appends the values physically stored on that node. For an internal node this walks
the descendant leaves and appends their stored values directly, avoiding the temporary vector that
`values(tree, node)` materializes.
"""
function appendvalues!(out, tree, node::Int)
    if iszero(firstchild(tree, node))
        append!(out, values(data(tree, node)))
    else
        for child in children(tree, node)
            appendvalues!(out, tree, child)
        end
    end
    return out
end

"""
    appendvalues!(out, tree, nodes)

Append all values contained in every node in `nodes` to `out`.
"""
function appendvalues!(out, tree, nodes)
    for node in nodes
        appendvalues!(out, tree, node)
    end
    return out
end

"""
    values(tree, node::Int)

Return the values contained in `node`'s subtree.

For a leaf, the stored value vector is returned directly. For an internal node,
the descendant leaf values are materialized into a new vector. Use
[`appendvalues!`](@ref), [`foreachvalue`](@ref), or [`anyvalue`](@ref) in
allocation-sensitive code.
"""
function values(tree, node::Int)
    iszero(firstchild(tree, node)) && return H2Trees.values(data(tree, node))
    return appendvalues!(Int[], tree, node)
end

"""
    values(tree, nodes)

Materialize all values contained in every node in `nodes`.
"""
function values(tree, nodes)
    return appendvalues!(Int[], tree, nodes)
end

"""
    values(data::Union{BoxData,BoundingBallData})

Return the value vector physically stored in a node-data payload.
"""
function values(data::Union{BoxData,BoundingBallData})
    return data.values
end

"""
    root(tree)

Return the root node id.
"""
function root(tree)
    return tree.root
end

"""
    level(tree, node)
    level(data)

Return the tree level of a node or node-data payload.
"""
function level(tree, nodeid::Int)
    return level(tree(nodeid).data)
end

function level(data::Union{BoxData,BoundingBallData})
    return data.level
end

"""
    nodes(tree)

Return the tree's node storage vector.
"""
function nodes(tree)
    return tree.nodes
end

"""
    lastnode(tree)

Return the largest global node id in `tree`.
"""
function lastnode(tree)
    return length(nodes(tree)) - root(tree) + 1
end

"""
    minimumlevel(tree)

Return the minimum level stored by `tree`.

This is usually the root level, but it need not be `1`.
"""
function minimumlevel(tree)
    return treeindex(tree).minlevel
end

"""
    levels(tree)

Return the level range covered by `tree`.
"""
function levels(tree)
    index = treeindex(tree)
    return index.minlevel:index.maxlevel
end

"""
    numberoflevels(tree)

Return the number of levels covered by `tree`.
"""
function numberoflevels(tree)
    return length(levels(tree))
end

"""
    parent(tree, node)
    nextsibling(tree, node)
    firstchild(tree, node)

Return topology links for `node`.
"""
function parent(tree, node::Int)
    return tree(node).parent
end

function nextsibling(tree, node::Int)
    return tree(node).nextsibling
end

function firstchild(tree, node::Int)
    return tree(node).firstchild
end

"""
    numberofnodes(tree)

Return the number of nodes stored by `tree`.
"""
function numberofnodes(tree)
    return length(tree.nodes)
end

"""
    isleaf(tree, node)

Return whether `node` has no children.
"""
function isleaf(tree, node::Int)
    return iszero(tree(node).firstchild)
end

"""
    findleafnode(tree, value::Int)

Return the leaf node containing `value`, or `0` when none is found.
"""
function findleafnode(tree, value::Int)
    for leaf in leaves(tree)
        (value ∈ H2Trees.values(data(tree, leaf))) && return leaf
    end
    return 0
end

"""
    children(tree, node=root(tree))

Return an iterator over the direct children of `node`.
"""
function children(tree, node::Int=root(tree))
    return ChildIterator(tree, node)
end

"""
    center(tree, node=root(tree))
    center(data)

Return a node or node-data center.
"""
function center(tree, nodeid::Int=root(tree))
    return center(tree(nodeid).data)
end

function center(data::Union{BoxData,BoundingBallData})
    return data.center
end

"""
    sector(tree, node)
    sector(data::BoxData)

Return the sector stored for a box-tree node.
"""
function sector(tree, node::Int)
    return sector(tree(node).data)
end

function sector(data::BoxData)
    return data.sector
end

"""
    halfsize(tree, node=root(tree))
    halfsize(data::BoxData)

Return the halfsize of a box-tree node.
"""
function halfsize(tree, nodeid::Int=root(tree))
    return halfsize(data(tree, nodeid))
end

function halfsize(data::BoxData)
    return data.halfsize
end

"""
    halfsizes(tree)

Return one representative halfsize per tree level.
"""
function halfsizes(tree)
    halfsizes = eltype(eltype(tree))[]
    for level in levels(tree)
        for node in nodesatlevel(tree, level)
            push!(halfsizes, halfsize(tree, node))
            break
        end
    end
    return halfsizes
end

"""
    minhalfsize(tree)

Return the halfsize at the deepest cached level.
"""
function minhalfsize(tree)
    nodesatlowestlevel = nodesatlevel(tree)[end]
    nodesatlowestlevel == Int[] && return 0

    return data(tree, nodesatlowestlevel[begin]).halfsize
end

"""
    radius(tree, node)
    radius(data::BoundingBallData)

Return the radius of a ball-tree node.
"""
function radius(tree, node::Int)
    return radius(tree(node).data)
end

function radius(data::BoundingBallData)
    return data.radius
end

"""
    treewithmorelevels(blocktree)

Return the side of a block tree with at least as many levels as the other side.
"""
function treewithmorelevels(tree)
    if length(levels(trialtree(tree))) >= length(levels(testtree(tree)))
        return trialtree(tree)
    else
        return testtree(tree)
    end
end

"""
    samelevelnodes(tree, node)

Return all nodes at `node`'s level.
"""
function samelevelnodes(tree, node::Int)
    return nodesatlevel(tree, level(tree, node))
end

"""
    nodesatlevel(tree)
    nodesatlevel(tree, level)

Return cached nodes grouped by level, or the nodes at one level.
"""
function nodesatlevel(tree)
    return treeindex(tree).nodes_by_level
end

function nodesatlevel(tree, level::Int)
    levels_ = levels(tree)
    level < levels_[begin] && return Int[]
    level > levels_[end] && return Int[]
    return nodesatlevel(tree)[leveltolevelid(tree, level)]
end

"""
    LevelIterator(tree, level::Int)

Return an iterator over nodes at `level`.
"""
function LevelIterator(tree, level::Int)
    return nodesatlevel(tree, level)
end

"""
    SameLevelIterator(tree, node::Int)

Return an iterator over nodes at the same level as `node`.
"""
function SameLevelIterator(tree, node::Int)
    return samelevelnodes(tree, node)
end

function _adjustnodesatlevels!(tree)
    return rebuildtreeindex!(tree)
end

"""
    rebuildtreeindex!(tree)

Recompute and replace the cached [`TreeIndex`](@ref) for `tree`.
"""
function rebuildtreeindex!(tree)
    # Atomic whole-value replacement: the tree struct is immutable and holds the index by
    # reference in a `Ref`, so this swaps the box's contents rather than mutating fields in place.
    tree.index[] = TreeIndex(tree)
    return tree
end

"""
    depthfirstnodes(tree)

Return the cached depth-first node order.
"""
function depthfirstnodes(tree)
    return treeindex(tree).depthfirstnodes
end

"""
    numberofvalues(tree, node::Int=root(tree))

Number of values contained in `node`'s subtree.

Walks the child subtrees directly, like [`appendvalues!`](@ref)/[`foreachvalue`](@ref)/
[`anyvalue`](@ref), rather than materializing a leaf or value vector via `leaves(tree, node)`/
`values(tree, node)`.
"""
function numberofvalues(tree, node::Int=root(tree))
    if iszero(firstchild(tree, node))
        return length(H2Trees.values(data(tree, node)))
    end
    nvals = 0
    for child in children(tree, node)
        nvals += numberofvalues(tree, child)
    end
    return nvals
end

"""
    valuesatnodes(tree; numberofvalues=H2Trees.numberofvalues(tree))

Return a vector mapping each stored value id to the leaves containing it.
"""
function valuesatnodes(tree; numberofvalues=numberofvalues(tree))
    nodes = Vector{Vector{Int}}(undef, numberofvalues)

    for leaf in leaves(tree)
        for value in H2Trees.values(data(tree, leaf))
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

"""
    nodesatvalues(tree, boxes=valuesatnodes(tree))

Group value ids by the set of leaf nodes that contain them.
"""
function nodesatvalues(tree, boxes=valuesatnodes(tree))
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

"""
    leveltolevelid(tree, level::Int)

Convert a tree level to the cached level-array index.
"""
function leveltolevelid(tree, level::Int)
    return level - minimumlevel(tree) + 1
end

"""
    levelindex(tree, node)

Return the cached level-array index for `node`.
"""
function levelindex(tree, node::Int)
    return leveltolevelid(tree, level(tree, node))
end

"""
    checkbalancedtree(tree)

Return whether all leaves sit at the same tree level.
"""
function checkbalancedtree(tree)
    leaflevel = level(tree, leaves(tree)[1])
    for node in leaves(tree)
        leaflevel != level(tree, node) && return false
    end
    return true
end

"""
    computevectorbuffers(tree, T)

Allocate one vector buffer per leaf, sized to that leaf's stored values.

For block trees, returns separate dictionaries for the test and trial sides.
"""
function computevectorbuffers(tree, T)
    return computevectorbuffers(tree, treetrait(tree), T)
end

function computevectorbuffers(tree, ::isBlockTree, ::Type{T}) where {T}
    return computevectorbuffers(testtree(tree), T), computevectorbuffers(trialtree(tree), T)
end

function computevectorbuffers(tree, ::AbstractTreeTrait, ::Type{T}) where {T}
    vectors = Vector{Vector{T}}(undef, length(leaves(tree)))

    for (i, leaf) in enumerate(leaves(tree))
        vectors[i] = Vector{T}(undef, length(H2Trees.values(data(tree, leaf))))
    end
    return Dict(zip(leaves(tree), vectors))
end

"""
    isuppertreelevel(tree, level)
    islowertreelevel(tree, level)
    isuppertreenode(tree, node)
    islowertreenode(tree, node)

Query whether a level or node belongs to the upper or lower part of a hybrid
tree.
"""
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
    return !isuppertreenode(tree, node)
end

include("trees/clustertrees.jl")
include("trees/TwoNTree.jl")
include("trees/treewrappers.jl")
include("trees/SimpleHybridTree.jl")
include("trees/BoundingBallTree.jl")
include("trees/BlockTree.jl")
include("trees/KMeansTree.jl")
include("trees/MetisTree.jl")
include("buildtree.jl")
include("showmethods.jl")

"""
    isgalerkinsymmetric(T)
    isgalerkinsymmetric(::Type{T})

Return whether a Galerkin operator type can be treated as symmetric.

Packages defining concrete operators can extend this hook. The fallback is
`false`.
"""
function isgalerkinsymmetric(T)
    return isgalerkinsymmetric(typeof(T))
end # requires BEAST to load

function isgalerkinsymmetric(::Type{T}) where {T}
    return false
end

export TwoNTree, BlockTree, SimpleHybridTree, KMeansTree

"""
    adjacencygraph(...)

Extension hook implemented by graph/visualization integrations.
"""
function adjacencygraph end # requires BEAST to load

# for backwards compatibility with julia versions below 1.9
if !isdefined(Base, :get_extension)
    include("../ext/H2BEASTTrees/H2BEASTTrees.jl")
    include("../ext/H2ParallelKMeansTrees/H2ParallelKMeansTrees.jl")
    include("../ext/H2PlotlyJSTrees/H2PlotlyJSTrees.jl")
end

end # module H2Trees
