
# Is point inside node #####################################################################

function isin(tree, node, point)
    return isin(tree, node, point, treetrait(tree))
end

# Distance measuring functions #############################################################

struct _IsNearFunctor{P}
    kwargs::P
end

function (f::_IsNearFunctor)(tree)
    return isnear(tree, treetrait(tree); f.kwargs...)
end

"""
    isnear
"""
function isnear(; kwargs...)
    return _IsNearFunctor(kwargs)
end

struct _IsNearNotBlockTreeFunctor{P}
    kwargs::P
end

function (f::_IsNearNotBlockTreeFunctor)(tree, testnode, trialnode)
    return isnear(tree, testnode, trialnode, treetrait(tree); f.kwargs...)
end

function isnear(tree, ::Any; kwargs...)
    return _IsNearNotBlockTreeFunctor(kwargs)
end
struct _IsFarNotBlockTreeFunctor{P}
    kwargs::P
end

function (f::_IsFarNotBlockTreeFunctor)(tree, testnode, trialnode)
    return !isnear(tree, testnode, trialnode, treetrait(tree); f.kwargs...)
end

function isfar(f::_IsNearNotBlockTreeFunctor)
    return _IsFarNotBlockTreeFunctor(f.kwargs)
end

struct _IsNearBlockTreeFunctor{P}
    kwargs::P
end

function (f::_IsNearBlockTreeFunctor)(testtree, trialtree, testnode, trialnode)
    return isnear(
        testtree,
        trialtree,
        testnode,
        trialnode,
        treetrait(testtree),
        treetrait(trialtree);
        f.kwargs...,
    )
end

struct _IsFarBlockTreeFunctor{P}
    kwargs::P
end

function (f::_IsFarBlockTreeFunctor)(testtree, trialtree, testnode, trialnode)
    return !isnear(
        testtree,
        trialtree,
        testnode,
        trialnode,
        treetrait(testtree),
        treetrait(trialtree);
        f.kwargs...,
    )
end

function isfar(f::_IsNearBlockTreeFunctor)
    return _IsFarBlockTreeFunctor(f.kwargs)
end

function isnear(tree, ::isBlockTree; kwargs...)
    return _IsNearBlockTreeFunctor(kwargs)
end

function isnear(tree, testnode::Int, trialnode::Int; kwargs...)
    return isnear(tree, testnode, trialnode, treetrait(tree), kwargs...)
end

function isnear(testtree, trialtree, testnode::Int, trialnode::Int)
    return isnear(
        testtree, trialtree, testnode, trialnode, treetrait(testtree), treetrait(trialtree)
    )
end

function isfar(tree, testnode::Int, trialnode::Int)
    return !isnear(tree, testnode::Int, trialnode::Int)
end

function isfar(testtree, trialtree, testnode::Int, trialnode::Int)
    return !isnear(
        testtree, trialtree, testnode, trialnode, treetrait(testtree), treetrait(trialtree)
    )
end

# TwoNTree #################################################################################

function isnear(
    tree,
    testnode::Int,
    trialnode::Int,
    ::isTwoNTree;
    additionalbufferboxes::Int=0,
    kwargs...,
)
    return isnearhalfsize(
        center(tree, testnode),
        center(tree, trialnode),
        halfsize(tree, testnode),
        additionalbufferboxes;
        kwargs...,
    )
end

function isnear(
    testtree,
    trialtree,
    testnode::Int,
    trialnode::Int,
    ::isTwoNTree,
    ::isTwoNTree;
    additionalbufferboxes::Int=0,
    kwargs...,
)
    return isnearhalfsize(
        center(testtree, testnode),
        center(trialtree, trialnode),
        max(halfsize(trialtree, trialnode), halfsize(testtree, testnode)),
        additionalbufferboxes;
        kwargs...,
    )
end

function isnearhalfsize(
    center_a::AbstractVector,
    center_b::AbstractVector,
    halfsize::T,
    additionalbufferboxes::Int;
    minhalfsizeadditionalbufferboxes::T=T(Inf),
    kwargs...,
) where {T}
    difference = center_a - center_b

    distancesquared = LinearAlgebra.dot(difference, difference)

    adbufferboxes = if minhalfsizeadditionalbufferboxes < halfsize
        0
    else
        additionalbufferboxes
    end

    return distancesquared <= (adbufferboxes + 1) * 12 * (1 + 100 * eps(T)) * halfsize^2
end

function isnearhalfsize(
    center_a::AbstractVector,
    center_b::AbstractVector,
    halfsize,
    additionalbufferboxes::Int,
    minhalfsize;
    kwargs...,
)
    halfsize < minhalfsize && return false
    return isnearhalfsize(center_a, center_b, halfsize, additionalbufferboxes; kwargs...)
end

function isin(tree, node, point, ::isTwoNTree)
    return isnearhalfsize(center(tree, node), point, halfsize(tree, node), 0, 0)
end

# BoundingBallTree #########################################################################

function isnear(tree, testnode::Int, trialnode::Int, ::isBoundingBallTree; kwargs...)
    return isnearradius(
        center(tree, testnode),
        center(tree, trialnode),
        radius(tree, testnode),
        radius(tree, trialnode);
        kwargs...,
    )
end

function isnear(
    testtree,
    trialtree,
    testnode::Int,
    trialnode::Int,
    ::isBoundingBallTree,
    ::isBoundingBallTree;
    kwargs...,
)
    return isnearradius(
        center(testtree, testnode),
        center(trialtree, trialnode),
        radius(testtree, testnode),
        radius(trialtree, trialnode),
        kwargs...,
    )
end

# η has to be >= 1
function isnearradius(
    center1::AbstractVector,
    center2::AbstractVector,
    radius1::T,
    radius2::T;
    η::T=one(T),
    kwargs...,
) where {T}
    @assert η >= one(T)
    difference = center1 - center2
    differencenorm = norm(difference)

    if (differencenorm + radius2 <= radius1)
        # ball2 is inside ball1
        return true

    elseif (differencenorm + radius1 <= radius2)
        # ball1 is inside ball2
        return true
    end

    maxradius = max(radius1, radius2)

    return differencenorm <= η * (1 + 10 * eps(T)) * (2 * maxradius)
end

function isin(tree, node, point, ::isBoundingBallTree)
    return isnearradius(
        center(tree, node), point, radius(tree, node), zero(radius(tree, node))
    )
end
