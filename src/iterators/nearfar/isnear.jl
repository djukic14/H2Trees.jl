
# Is point inside node #####################################################################

"""
    isin(tree, node, point)

Return whether `point` lies in `node` for `tree`.
"""
function isin(tree, node, point)
    return isin(tree, node, point, treetrait(tree))
end

# Distance measuring functions #############################################################

struct IsNearFunctor{P}
    kwargs::P
end

function (f::IsNearFunctor)(tree)
    return isnear(tree, treetrait(tree); f.kwargs...)
end

function isnear(; kwargs...)
    return IsNearFunctor(kwargs)
end

struct IsNearNotBlockTreeFunctor{P}
    kwargs::P
end

function (f::IsNearNotBlockTreeFunctor)(tree, testnode, trialnode)
    return isnear(tree, testnode, trialnode; f.kwargs...)
end

function isnear(tree, ::Any; kwargs...)
    return IsNearNotBlockTreeFunctor(kwargs)
end

struct IsFarNotBlockTreeFunctor{P}
    kwargs::P
end

function (f::IsFarNotBlockTreeFunctor)(tree, testnode, trialnode)
    return !isnear(tree, testnode, trialnode, treetrait(tree); f.kwargs...)
end

function isfar(f::IsNearNotBlockTreeFunctor)
    return IsFarNotBlockTreeFunctor(f.kwargs)
end
struct IsNearBlockTreeFunctor{P}
    kwargs::P
end

function (f::IsNearBlockTreeFunctor)(testtree, trialtree, testnode, trialnode)
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

function isfar(f::IsNearBlockTreeFunctor)
    return _IsFarBlockTreeFunctor(f.kwargs)
end

function isnear(tree, ::isBlockTree; kwargs...)
    return IsNearBlockTreeFunctor(kwargs)
end

"""
    isnear(tree, testnode::Int, trialnode::Int; minlevel=level(tree, root(tree)), kwargs...)

Return whether two nodes in a single tree are near.

If `level(tree, testnode) < minlevel`, this returns `true` immediately;
otherwise it forwards to trait-dispatched `isnear`.
"""
function isnear(
    tree, testnode::Int, trialnode::Int; minlevel::Int=level(tree, root(tree)), kwargs...
)
    level(tree, testnode) < minlevel && return true
    return isnear(tree, testnode, trialnode, treetrait(tree); kwargs...)
end

"""
    isnear(testtree, trialtree, testnode::Int, trialnode::Int; minlevel=level(testtree, root(testtree)), kwargs...)

Return whether two nodes, one from each tree, are near.

If either node level is below `minlevel`, this returns `true` immediately;
otherwise it forwards to trait-dispatched `isnear`.
"""
function isnear(
    testtree,
    trialtree,
    testnode::Int,
    trialnode::Int;
    minlevel::Int=level(testtree, root(testtree)),
    kwargs...,
)
    H2Trees.level(testtree, testnode) < minlevel && return true
    H2Trees.level(trialtree, trialnode) < minlevel && return true

    return isnear(
        testtree,
        trialtree,
        testnode,
        trialnode,
        treetrait(testtree),
        treetrait(trialtree);
        kwargs...,
    )
end

"""
    isfar(tree, testnode::Int, trialnode::Int)

Return the logical negation of `isnear(tree, testnode, trialnode)`.
"""
function isfar(tree, testnode::Int, trialnode::Int)
    return !isnear(tree, testnode::Int, trialnode::Int)
end

"""
    isfar(testtree, trialtree, testnode::Int, trialnode::Int)

Return the logical negation of the trait-dispatched two-tree `isnear` predicate.
"""
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
    kwargs...,
) where {T}
    difference = center_a - center_b

    distancesquared = LinearAlgebra.dot(difference, difference)

    return distancesquared <=
           (additionalbufferboxes + 1) * 12 * (1 + 100 * eps(T)) * halfsize^2
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
