
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

"""
    isnear(; kwargs...)

Return a near-predicate factory.

The returned object resolves to the appropriate one-tree or block-tree predicate
when called with a tree. Keyword arguments are forwarded to the concrete
geometric predicate.
"""
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
    level(testtree, testnode) < minlevel && return true
    level(trialtree, trialnode) < minlevel && return true

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
    return isneargap(
        center(tree, testnode),
        center(tree, trialnode),
        halfsize(tree, testnode),
        halfsize(tree, trialnode),
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
    # Keep each side's own halfsize; independently built trees can carry
    # different box sizes at the same raw level.
    return isneargap(
        center(testtree, testnode),
        center(trialtree, trialnode),
        halfsize(testtree, testnode),
        halfsize(trialtree, trialnode),
        additionalbufferboxes;
        kwargs...,
    )
end

"""
    NEARGAPRELTOL

Relative tolerance for the box-gap near predicate.

This is only the numerical slack around the admissibility margin. The margin
itself is [`DEFAULTNEARGAPBOXES`](@ref).
"""
const NEARGAPRELTOL = 1e-3

"""
    DEFAULTNEARGAPBOXES

Default admissibility margin for [`isneargap`](@ref), measured in units of the
smaller box's halfsize.

The value `1.0` is the midpoint between touching boxes and the next
non-neighbour on a common grid: those gaps are `0` and `2h`, so every threshold
in `(0, 2h)` preserves the usual 3^d near neighbourhood. The midpoint also gives
independently built block-tree sides a real separation buffer without placing
the boundary on a common-grid lattice value. `additionalbufferboxes` adds more
half-sizes on top.
"""
const DEFAULTNEARGAPBOXES = 1.0

"""
    isneargap(center_a, center_b, halfsize_a, halfsize_b, additionalbufferboxes;
              gapboxes=DEFAULTNEARGAPBOXES, reltol=NEARGAPRELTOL)

Box-gap near predicate: near when the axis-aligned gap between the two boxes is at most
`(gapboxes + additionalbufferboxes) * min(halfsize)`, plus a small relative tolerance.

Preferred over [`isnearhalfsize`](@ref) for admissibility. It measures the actual box gap using
each side's own halfsize, so independently built trees need not share a grid. See
[`DEFAULTNEARGAPBOXES`](@ref) for why the default margin is one half-size.
"""
function isneargap(
    center_a::AbstractVector,
    center_b::AbstractVector,
    halfsize_a::T,
    halfsize_b,
    additionalbufferboxes::Int;
    gapboxes=DEFAULTNEARGAPBOXES,
    reltol=NEARGAPRELTOL,
    kwargs...,
) where {T}
    reach = halfsize_a + halfsize_b
    gapsquared = zero(T)
    for i in eachindex(center_a)
        gap = abs(center_a[i] - center_b[i]) - reach
        gap > 0 && (gapsquared += gap * gap)
    end
    # The margin follows the smaller box; the tolerance follows the larger box
    # so it stays a true relative slack.
    allowance = (gapboxes + additionalbufferboxes) * min(halfsize_a, halfsize_b)
    return sqrt(gapsquared) <= allowance + reltol * max(halfsize_a, halfsize_b)
end

"""
    isnearhalfsize(center_a, center_b, halfsize, additionalbufferboxes)

LEGACY centre-distance near predicate: near when the centres lie within `sqrt((n+1)*12)*halfsize`.

Retained as a named predicate, but not used for default admissibility. It
assumes both boxes share one `halfsize`; use [`isneargap`](@ref) for
independently built trees.
"""
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
    nodecenter = center(tree, node)
    nodehalfsize = halfsize(tree, node)
    tolerance = 100 * eps(typeof(nodehalfsize)) * nodehalfsize

    for i in eachindex(nodecenter)
        abs(point[i] - nodecenter[i]) <= nodehalfsize + tolerance || return false
    end

    return true
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
        radius(trialtree, trialnode);
        kwargs...,
    )
end

"""
    isnearradius(center1, center2, radius1, radius2; η=1)

Ball near predicate.

Two balls are near when one contains the other, or when their center distance is
at most `η * (radius1 + radius2)` up to numerical slack. `η` must be at least one.
"""
function isnearradius(
    center1::AbstractVector,
    center2::AbstractVector,
    radius1::T,
    radius2::T;
    η::T=one(T),
    kwargs...,
) where {T}
    # Two bounding balls are NEAR when the distance between their centres does not exceed
    # `η * (radius1 + radius2)`; equivalently, when their *gap* `‖c1 - c2‖ - r1 - r2` is at most
    # `(η - 1) * (r1 + r2)`. This sum-of-radii form is used instead of `2 * max(radius)` because
    # it makes the near/far classification MONOTONE across tree levels, which the well-separated
    # translation partition (and `testwellseparatedness`) requires:
    #
    #   * a parent ball encloses each of its children, so `r_parent >= r_child` on both sides
    #     and the gap can only shrink going up the tree (triangle inequality);
    #   * hence if two child nodes are near (small gap), their parents are necessarily near too.
    #
    # The old `2 * max(radius)` form has an effective margin of `|r1 - r2|`, which is NOT
    # monotone up the tree, so with minimal SEBB parent radii a near child pair could sit under
    # a far parent pair: a double-counted translation. The two forms coincide for equal radii
    # (the common same-level case), so this preserves the intended behaviour while fixing the
    # monotonicity defect at the admissibility boundary.
    #
    # Note the argument above needs only `r_parent >= r_child` (CONTAINMENT) and not that the
    # parent ball is minimal. So it survives SEBB's approximate fallback unchanged: a fallback
    # ball's radius is the measured enclosing radius at its centre, so it still contains its
    # children. An approximate parent is merely LARGER, which can only widen the near set, never
    # break monotonicity.
    if η < one(T)
        throw(ArgumentError("η must be greater than or equal to 1, got $η"))
    end
    difference = center1 - center2
    differencenorm = norm(difference)

    if (differencenorm + radius2 <= radius1)
        # ball2 is inside ball1
        return true

    elseif (differencenorm + radius1 <= radius2)
        # ball1 is inside ball2
        return true
    end

    return differencenorm <= η * (1 + 100 * eps(T)) * (radius1 + radius2)
end

function isin(tree, node, point, ::isBoundingBallTree)
    return isnearradius(
        center(tree, node), point, radius(tree, node), zero(radius(tree, node))
    )
end
