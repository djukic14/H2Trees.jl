"""
    AbstractTreeTrait

Supertype for tree-dispatch tags.

Algorithms that depend on tree geometry dispatch on these tags instead of on
concrete tree types. New tree types can participate by defining
[`treetrait`](@ref).
"""
abstract type AbstractTreeTrait end

"""
    isTwoNTree <: AbstractTreeTrait

Trait tag for box-shaped trees, such as [`TwoNTree`](@ref).
"""
struct isTwoNTree <: AbstractTreeTrait end

"""
    isBoundingBallTree <: AbstractTreeTrait

Trait tag for ball-shaped trees, such as [`BoundingBallTree`](@ref) and
[`KMeansTree`](@ref).
"""
struct isBoundingBallTree <: AbstractTreeTrait end

"""
    isBlockTree <: AbstractTreeTrait

Trait tag for a [`BlockTree`](@ref), i.e. a test/trial tree pair.
"""
struct isBlockTree <: AbstractTreeTrait end

"""
    isAnyTree <: AbstractTreeTrait

Fallback trait used when code should exercise the most general tree path.
"""
struct isAnyTree <: AbstractTreeTrait end

"""
    treetrait(tree)

Return the [`AbstractTreeTrait`](@ref) tag for `tree`.
"""
function treetrait(tree)
    return treetrait(typeof(tree))
end
