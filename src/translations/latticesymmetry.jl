"""
    LatticeSymmetry{D}

One symmetry of the cubic lattice: a permutation of the axes, each optionally flipped in sign.

`permutation[i]` says which axis the `i`-th output component comes from, and `signs[i]` whether it
is negated:

    (S q)[i] = signs[i] * q[permutation[i]]

That is [`applysymmetry`](@ref). There are `D! * 2^D` such symmetries: 8 in 2-D, 48 in 3-D.

These maps preserve dot products, which is what lets one translation serve several offsets. If a
translation operator depends on direction only through `k̂ ⋅ d`, then `(S k̂) ⋅ (S d) == k̂ ⋅ d`,
and so `T(k̂, S d) == T(S⁻¹ k̂, d)`.
"""
struct LatticeSymmetry{D}
    permutation::NTuple{D,Int8}
    signs::NTuple{D,Int8}

    function LatticeSymmetry(permutation::NTuple{D,Int8}, signs::NTuple{D,Int8}) where {D}
        for p in permutation
            (1 <= p <= D) ||
                error("LatticeSymmetry: permutation entry $(p) is outside 1:$(D).")
        end
        length(unique(permutation)) == D ||
            error("LatticeSymmetry: $(permutation) is not a permutation of 1:$(D).")
        all(s -> s == 1 || s == -1, signs) ||
            error("LatticeSymmetry: signs $(signs) must all be +1 or -1.")
        return new{D}(permutation, signs)
    end
end

function LatticeSymmetry(
    permutation::NTuple{D,<:Integer}, signs::NTuple{D,<:Integer}
) where {D}
    return LatticeSymmetry(Int8.(permutation), Int8.(signs))
end

function Base.:(==)(a::LatticeSymmetry{D}, b::LatticeSymmetry{D}) where {D}
    return a.permutation == b.permutation && a.signs == b.signs
end

function Base.hash(s::LatticeSymmetry, h::UInt)
    return hash(s.signs, hash(s.permutation, hash(:LatticeSymmetry, h)))
end

"""
    identitysymmetry(Val(D))

The identity element of the `D`-dimensional lattice symmetry group.
"""
function identitysymmetry(::Val{D}) where {D}
    return LatticeSymmetry(ntuple(i -> Int8(i), Val(D)), ntuple(_ -> Int8(1), Val(D)))
end

# `Number`, not `Integer`: a lattice symmetry is an orthogonal map, and consumers apply it to
# direction vectors as well as to integer offsets. MLFMA derives its angular permutations by
# asking where a symmetry sends each sample point. The element type is deliberately not a type
# parameter (`NTuple{0,T}` would leave it unbound, which Aqua flags); `Int8 * Number` promotes to
# the argument's own type, which is what is wanted anyway.
"""
    applysymmetry(symmetry, q)

Apply `symmetry` to `q`: component `i` of the result is `signs[i] * q[permutation[i]]`.

Takes and returns a tuple. Use [`applysymmetry!`](@ref) if you are working with vectors, as the
translation build loop does.
"""
function applysymmetry(symmetry::LatticeSymmetry{D}, q::NTuple{D,<:Number}) where {D}
    return ntuple(i -> symmetry.signs[i] * q[symmetry.permutation[i]], Val(D))
end

"""
    applysymmetry!(destination, symmetry, q)

The vector form of [`applysymmetry`](@ref), writing the result into `destination`.

`destination` and `q` must not be the same array.
"""
function applysymmetry!(
    destination::AbstractVector, symmetry::LatticeSymmetry{D}, q::AbstractVector
) where {D}
    for i in 1:D
        destination[i] = symmetry.signs[i] * q[symmetry.permutation[i]]
    end
    return destination
end

"""
    composesymmetry(a, b)

The symmetry that applies `b` first and then `a`:

    applysymmetry(composesymmetry(a, b), q) == applysymmetry(a, applysymmetry(b, q))
"""
function composesymmetry(a::LatticeSymmetry{D}, b::LatticeSymmetry{D}) where {D}
    permutation = ntuple(i -> b.permutation[a.permutation[i]], Val(D))
    signs = ntuple(i -> a.signs[i] * b.signs[a.permutation[i]], Val(D))
    return LatticeSymmetry(permutation, signs)
end

"""
    determinant(symmetry) -> Int

The determinant of `symmetry`: `+1` for a rotation, `-1` for a reflection.

Computed exactly, as the permutation's parity times the product of the signs.

You need this when a symmetry acts on field components rather than just on positions. A
pseudo-tensor (anything built from a curl or a cross product, such as a double-layer kernel)
picks up the determinant where an ordinary tensor does not:

    tensor         G(Sx) == S G(x) Sᵀ
    pseudo-tensor  K(Sx) == det(S) · S K(x) Sᵀ

Both are 3x3, so their shape does not tell you which you have. Using the wrong rule flips the sign
of the result on every reflection and is correct on every rotation. Reflections are half of any
lattice symmetry group, so a test that only tries the identity or a rotation will not catch it.
"""
function determinant(symmetry::LatticeSymmetry{D}) where {D}
    # Permutation parity by counting inversions. D is 2 or 3, so this is cheaper than anything
    # cleverer and needs no allocation.
    inversions = 0
    for i in 1:D, j in (i + 1):D
        symmetry.permutation[i] > symmetry.permutation[j] && (inversions += 1)
    end
    parity = iseven(inversions) ? 1 : -1
    return parity * prod(Int, symmetry.signs)
end

"""
    inversesymmetry(symmetry)

The inverse of `symmetry`, satisfying `composesymmetry(inversesymmetry(s), s) == identitysymmetry(Val(D))`.
"""
function inversesymmetry(symmetry::LatticeSymmetry{D}) where {D}
    inversepermutation = ntuple(Val(D)) do i
        return Int8(findfirst(==(Int8(i)), symmetry.permutation))
    end
    signs = ntuple(i -> symmetry.signs[inversepermutation[i]], Val(D))
    return LatticeSymmetry(inversepermutation, signs)
end

"""
    AbstractLatticeSymmetryGroup

Supertype for the symmetry policies accepted by [`symmetrygroup`](@ref).
"""
abstract type AbstractLatticeSymmetryGroup end

"""
    NoSymmetry <: AbstractLatticeSymmetryGroup

The one-element group. Every offset is its own representative, so this gives the same result as
plain direction deduplication.
"""
struct NoSymmetry <: AbstractLatticeSymmetryGroup end

"""
    OppositeSymmetry <: AbstractLatticeSymmetryGroup

The two-element group `{identity, -identity}`, treating `q` and `-q` as the same offset.

The weakest useful choice, and the safest: it asks only that the angular sampling be closed under
`k̂ -> -k̂`.
"""
struct OppositeSymmetry <: AbstractLatticeSymmetryGroup end

"""
    AxisPreservingSymmetry{A} <: AbstractLatticeSymmetryGroup

The symmetries that keep axis `A` on itself, possibly flipped, and let the other axes permute and
change sign freely. `(D-1)! * 2^D` elements: 16 of the 48 in 3-D when `A = 3`.
"""
struct AxisPreservingSymmetry{A} <: AbstractLatticeSymmetryGroup end

AxisPreservingSymmetry(axis::Int) = AxisPreservingSymmetry{axis}()

"""
    FullLatticeSymmetry <: AbstractLatticeSymmetryGroup

All `D! * 2^D` symmetries.

Only use this if the representation is closed under every rotation and reflection of the cube.
"""
struct FullLatticeSymmetry <: AbstractLatticeSymmetryGroup end

# All permutations of 1:n in lexicographic order, so element ordering (and hence every
# symmetry ID) is reproducible across runs and machines.
function _permutations(n::Int)
    n == 0 && return [Int8[]]
    out = Vector{Int8}[]
    for first in 1:n, rest in _permutations(n - 1)
        push!(out, Int8[first; Int8[r >= first ? r + 1 : r for r in rest]])
    end
    return out
end

function _allsymmetries(::Val{D}) where {D}
    out = LatticeSymmetry{D}[]
    for permutation in _permutations(D),
        signs in Iterators.product(ntuple(_ -> (1, -1), Val(D))...)

        push!(out, LatticeSymmetry(NTuple{D,Int}(permutation), signs))
    end
    return out
end

_inpolicy(::FullLatticeSymmetry, ::LatticeSymmetry) = true
_inpolicy(::NoSymmetry, s::LatticeSymmetry{D}) where {D} = s == identitysymmetry(Val(D))
function _inpolicy(::OppositeSymmetry, s::LatticeSymmetry{D}) where {D}
    return s.permutation == identitysymmetry(Val(D)).permutation &&
           (all(==(Int8(1)), s.signs) || all(==(Int8(-1)), s.signs))
end
function _inpolicy(::AxisPreservingSymmetry{A}, s::LatticeSymmetry{D}) where {A,D}
    1 <= A <= D || error(
        "AxisPreservingSymmetry{$(A)} is not usable in $(D) dimensions; the axis must lie in 1:$(D).",
    )
    return s.permutation[A] == Int8(A)
end

"""
    SymmetryGroup{D,P}

A policy's symmetries in `D` dimensions, plus the lookup table canonicalization needs.

`elements[1]` is the identity and the order is reproducible, so an index into `elements` is a
stable symmetry ID, small enough to store once per interaction (a `UInt8` is always enough).
`inverseindex[k]` is where `inversesymmetry(elements[k])` sits. It always exists, because a policy
names a subgroup and subgroups are closed under inversion.

Build one with [`symmetrygroup`](@ref).
"""
struct SymmetryGroup{D,P<:AbstractLatticeSymmetryGroup}
    policy::P
    elements::Vector{LatticeSymmetry{D}}
    inverseindex::Vector{Int}
end

Base.length(group::SymmetryGroup) = length(group.elements)
# `Integer`, not `Int`: consumers store symmetry IDs in the narrowest type that fits (`UInt8`
# suffices for every group in every dimension), and should not have to widen to index the group.
Base.getindex(group::SymmetryGroup, i::Integer) = group.elements[i]
Base.eachindex(group::SymmetryGroup) = eachindex(group.elements)
Base.iterate(group::SymmetryGroup, state...) = iterate(group.elements, state...)
Base.eltype(::Type{SymmetryGroup{D,P}}) where {D,P} = LatticeSymmetry{D}

"""
    symmetrygroup(policy, Val(D))

Build `policy`'s symmetries in `D` dimensions as a [`SymmetryGroup`](@ref).

Checks that the policy really is a subgroup: identity first, closed under composition and
inversion. A policy that is not closed would not fail outright. It would quietly break
canonicalization instead, because two offsets in the same orbit could reduce to different
representatives.

```jldoctest
julia> length(H2Trees.symmetrygroup(H2Trees.FullLatticeSymmetry(), Val(3)))
48

julia> length(H2Trees.symmetrygroup(H2Trees.AxisPreservingSymmetry(3), Val(3)))
16
```
"""
function symmetrygroup(policy::AbstractLatticeSymmetryGroup, ::Val{D}) where {D}
    elements = filter(s -> _inpolicy(policy, s), _allsymmetries(Val(D)))

    elements[1] == identitysymmetry(Val(D)) || error(
        "symmetrygroup: $(policy) does not contain the identity as its first element."
    )

    index = Dict(s => k for (k, s) in enumerate(elements))
    inverseindex = Vector{Int}(undef, length(elements))
    for (k, s) in enumerate(elements)
        inverse = inversesymmetry(s)
        haskey(index, inverse) || error(
            "symmetrygroup: $(policy) is not closed under inversion in $(D) dimensions."
        )
        inverseindex[k] = index[inverse]
    end
    for a in elements, b in elements
        haskey(index, composesymmetry(a, b)) || error(
            "symmetrygroup: $(policy) is not closed under composition in $(D) dimensions.",
        )
    end

    return SymmetryGroup{D,typeof(policy)}(policy, elements, inverseindex)
end

# Lexicographic comparison, used to pick one representative per orbit. Any total order would do;
# this one puts the full group's representatives in the familiar q1 >= q2 >= ... >= qD >= 0
# region, which makes canonical offsets readable in diagnostics.
function _islarger(a::NTuple{D}, b::NTuple{D}) where {D}
    for i in 1:D
        a[i] == b[i] || return a[i] > b[i]
    end
    return false
end

function _islarger(a::AbstractVector, b::AbstractVector, indices)
    for i in indices
        a[i] == b[i] || return a[i] > b[i]
    end
    return false
end

"""
    canonicalizetranslation(q, group) -> (canonical, symmetryid)

Reduce `q` to the representative of its orbit under `group`, and return the ID of the symmetry
that maps that representative back to `q`:

    q == applysymmetry(group[symmetryid], canonical)

The stored symmetry goes from canonical to actual, not the other way round. That is the direction
a consumer wants: it holds the one translation built for `canonical` and needs the one for `q`. If
the operator depends on direction through `k̂ ⋅ d`, the consumer applies the inverse to its
directions:

    T(k̂, q) == T(applysymmetry(inversesymmetry(group[symmetryid]), k̂), canonical)

Worth knowing when testing this: the antipodal map is its own inverse, so `OppositeSymmetry` alone
cannot tell the two directions apart, and a reversed convention would still pass.

Everything in one orbit reduces to the same representative, which is what makes the result usable
as a deduplication key:

    canonicalizetranslation(applysymmetry(R, q), group)[1] ==
        canonicalizetranslation(q, group)[1]        for every R in group

Ties go to the lowest-numbered symmetry, so offsets with small orbits (those lying on a symmetry
axis or plane, like `(3, 0, 0)` or `(3, 3, 3)`) still give a deterministic answer.

```jldoctest
julia> group = H2Trees.symmetrygroup(H2Trees.FullLatticeSymmetry(), Val(3));

julia> canonical, id = H2Trees.canonicalizetranslation((-2, 3, -1), group);

julia> canonical
(3, 2, 1)

julia> H2Trees.applysymmetry(group[id], canonical)
(-2, 3, -1)
```
"""
function canonicalizetranslation(q::NTuple{D,<:Integer}, group::SymmetryGroup{D}) where {D}
    canonical = q
    best = 1
    for k in 2:length(group)
        candidate = applysymmetry(group[k], q)
        if _islarger(candidate, canonical)
            canonical = candidate
            best = k
        end
    end
    return canonical, group.inverseindex[best]
end

"""
    canonicalizetranslation!(canonical, scratch, q, group) -> symmetryid

The vector form of [`canonicalizetranslation`](@ref): writes the representative into `canonical`
and uses `scratch` as working space. Both need at least `D` entries, and neither may be the same
array as `q`.

This is what the translation build loop uses, since it keeps its offsets in reused vectors rather
than tuples.
"""
function canonicalizetranslation!(
    canonical::AbstractVector,
    scratch::AbstractVector,
    q::AbstractVector,
    group::SymmetryGroup{D},
) where {D}
    indices = 1:D
    copyto!(view(canonical, indices), view(q, indices))
    best = 1
    for k in 2:length(group)
        applysymmetry!(scratch, group[k], q)
        if _islarger(scratch, canonical, indices)
            copyto!(view(canonical, indices), view(scratch, indices))
            best = k
        end
    end
    return group.inverseindex[best]
end

"""
    symmetryorbit(q, group)

The set of distinct offsets reachable from `q` under `group`.

Orbits are smaller than the group for offsets fixed by some of its elements, which is why the
achievable translation reduction is measured rather than read off the group order.
"""
function symmetryorbit(q::NTuple{D,<:Integer}, group::SymmetryGroup{D}) where {D}
    return Set(applysymmetry(s, q) for s in group)
end
