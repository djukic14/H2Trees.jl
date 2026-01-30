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
