abstract type AbstractTreeTrait end

# something like an Octree, Quadtree
struct isTwoNTree <: AbstractTreeTrait end

struct isBoundingBallTree <: AbstractTreeTrait end

struct isBlockTree <: AbstractTreeTrait end

# for a tree not fitting in any other category (or the classification is avoided)
# mainly used for testing where the most general method version is tested
struct isAnyTree <: AbstractTreeTrait end

function treetrait(tree)
    return treetrait(typeof(tree))
end

abstract type UniquePointsTrait end

struct UniquePoints <: UniquePointsTrait end

# in case that some points are allowed to be in multiple nodes
struct NonUniquePoints <: UniquePointsTrait end

function uniquepointstreetrait(::Any)
    return UniquePoints()
end

function arepointsunique(tree)
    return arepointsunique(tree, treetrait(tree))
end

function arepointsunique(tree, ::AbstractTreeTrait)
    return uniquepointstreetrait(tree) isa UniquePoints
end

function arepointsunique(tree, ::isBlockTree)
    return arepointsunique(testtree(tree)) && arepointsunique(trialtree(tree))
end
