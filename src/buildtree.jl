# Canonical builder-only construction API. Public methods accept a `builder`
# keyword and forward to positional helpers so dispatch can use the builder type.

"""
    buildtree(input; builder=TwoNTreeBuilder())

Build a tree from one input object.

The builder selects the tree type and carries all construction settings.
Supported builders include [`TwoNTreeBuilder`](@ref),
[`BoundingBallTreeBuilder`](@ref), [`KMeansTreeBuilder`](@ref), and
[`SimpleHybridTreeBuilder`](@ref). The SimpleHybrid builder expects an existing
tree as input; the others expect a point collection.
"""
function buildtree(points; builder=TwoNTreeBuilder())
    return _buildtree(points, builder)
end

"""
    buildtree(testpoints, trialpoints; builder=BlockTreeBuilder())

Build a [`BlockTree`](@ref) from test and trial point collections.

The `BlockTreeBuilder` carries the two side-specific [`TwoNTreeBuilder`](@ref)s.
"""
function buildtree(testpoints, trialpoints; builder=BlockTreeBuilder())
    return _buildtree(testpoints, trialpoints, builder)
end

"""
    buildtree(points, graph, weights; builder=MetisTreeBuilder())

Build a [`MetisTree`](@ref) from points and a partitioning graph.

`weights` are the vertex weights used by METIS. The builder carries the
partitioning and tree-construction settings.
"""
function buildtree(points, graph, weights; builder=MetisTreeBuilder())
    return _buildtree(points, graph, weights, builder)
end

"""
    buildforest(points, graph, weights; builder=MetisForestBuilder())

Construct a [`MetisForest`](@ref), one tree per connected component of `graph`.
"""
function buildforest(points, graph, weights; builder=MetisForestBuilder())
    return _buildforest(points, graph, weights, builder)
end

_buildtree(points, builder::TwoNTreeBuilder) = TwoNTree(points; builder=builder)
function _buildtree(points, builder::BoundingBallTreeBuilder)
    return BoundingBallTree(points; builder=builder)
end
_buildtree(points, builder::KMeansTreeBuilder) = KMeansTree(points; builder=builder)
function _buildtree(tree, builder::SimpleHybridTreeBuilder)
    return SimpleHybridTree(tree; builder=builder)
end

function _buildtree(testpoints, trialpoints, builder::BlockTreeBuilder)
    return BlockTree(testpoints, trialpoints; builder=builder)
end

function _buildtree(points, graph, weights, builder::MetisTreeBuilder)
    return MetisTree(points, graph, weights; builder=builder)
end

function _buildforest(points, graph, weights, builder::MetisForestBuilder)
    return MetisForest(points, graph, weights; builder=builder)
end

function _buildtree(points, builder)
    return throw(
        ArgumentError(
            "buildtree does not know how to build from a single positional argument with a " *
            "$(typeof(builder)); pass a matching builder or use the multi-argument form",
        ),
    )
end
