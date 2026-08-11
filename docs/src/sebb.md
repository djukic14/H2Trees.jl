# Smallest Enclosing Ball of Balls (`SEBB`)

`H2Trees.SEBB` is a self-contained submodule that solves the *smallest enclosing ball of balls*
(SEBB) problem: given a finite set of input balls, find the unique ball of minimal radius that
contains all of them. It solves this exactly whenever the configuration is resolvable at the
caller's precision, and otherwise returns a strictly enclosing approximation rather than
failing (see [Degenerate inputs never raise](@ref)). The submodule is intentionally isolated.

## Problem

For input balls ``B_i = B(c_i, r_i)`` with ``c_i \in \mathbb{R}^d`` and ``r_i \ge 0``, SEBB
finds the ball ``B(c^\star, R^\star)`` with

```math
R^\star = \min_{c \in \mathbb{R}^d} \max_i \left( \lVert c - c_i \rVert + r_i \right).
```

A candidate ball ``B(c, R)`` contains input ``B_i`` exactly when
``\lVert c - c_i \rVert + r_i \le R``.

## What the solution looks like

The figure below shows one enclosing ball per supported dimension, with **every input ball
touching the enclosing ball from the inside**. Use the buttons to switch between the 1D, 2D and
3D cases:

```@raw html
<object data="../assets/plots/sebb_tangent.html" type="text/html"  style="width:100%; height:55vh;"> </object>
```

## Supported inputs

- ambient dimensions ``d = 1, 2, 3`` (other dimensions are rejected);
- `Float32`, `Float64`, and `BigFloat` (any `AbstractFloat`);
- finite centers;
- finite, **nonnegative** radii: positive radii are the primary use case, and zero radii are
  accepted as degenerate point balls (negative radii are rejected).

## Constructing balls

```@example sebb
using H2Trees
using StaticArrays
const SEBB = H2Trees.SEBB

b = SEBB.Ball(SVector(1.0, 2.0), 0.5)
(SEBB.center(b), SEBB.radius(b))
```

Integer or mixed scalar input is promoted to a common floating-point type
(`SEBB.Ball([1, 2, 3], 4)` builds a `Ball{3,Float64}`).

## Direct use (no tree required)

```@example sebb
balls = [
    SEBB.Ball(SVector(0.0, 0.0), 1.0),
    SEBB.Ball(SVector(3.0, 0.0), 0.5),
    SEBB.Ball(SVector(1.0, 2.0), 0.8),
]
ball = SEBB.smallest_enclosing_ball(balls)
(SEBB.center(ball), SEBB.radius(ball))
```

A convenience method accepts parallel `centers` and `radii` collections:

```@example sebb
centers = [SVector(0.0, 0.0), SVector(2.0, 0.0)]
radii = [1.0, 1.0]
SEBB.smallest_enclosing_ball(centers, radii)
```

## Algorithm

For dimension `d ∈ {1,2,3}`, the optimum is determined by at most `d + 1`
support balls with affinely independent centers.

Enumerate all supports of size `1:min(d + 1, n)`:

- size 1: return the ball itself;
- size 2: use the closed-form two-ball solution;
- size `m ≥ 3`: solve an `(m - 1) × (m - 1)` Gram system, then a quadratic
  equation for the enclosing radius.

Accept a candidate only if its center lies in the convex hull of the support
centers, the support balls are tangent within tolerance, and it encloses all
input balls.

Return the feasible candidate with smallest radius, using a deterministic
tie-breaker.

The cost is `O(n^(d+1))` for fixed `d`, which is appropriate for small sets such
as the child balls of a tree node.

## Numerical tolerance caveat

All comparisons route through one centralized, scale-aware tolerance policy that uses local
geometric scales (radii and relative distances), never absolute coordinate magnitudes, so the
result is translation invariant. The returned floating-point ball is inflated by any tiny
positive containment residual so that it is guaranteed to enclose every input under the same
arithmetic used by the caller. Tolerances are configurable through the `atol`/`rtol` keyword
arguments.

## Degenerate inputs never raise

`smallest_enclosing_ball` runs deep inside tree construction, where an exception over a
rounding accident would destroy an entire build. It therefore never fails on a valid, nonempty
input.

If no support set produces a numerically valid tangent ball at the caller's precision, the
solver emits a single warning and returns an *approximate* ball instead. That ball is still
guaranteed to contain every input (its radius is the measured enclosing radius at its centre,
so containment holds by construction rather than by tolerance) but it is not certified
minimal.

The public `smallest_enclosing_ball` API returns only the enclosing ball. Tests and internal
diagnostics can call `_smallest_enclosing_ball_with_certificate` to distinguish the exact and
fallback paths: an approximate result has an empty support, so `isexactresult(result) == false`.
For a bounding-ball tree the distinction is harmless (a node ball only has to *contain* its
children) but it is worth knowing when the exact radius matters, in which case a wider float
type usually resolves the configuration.

## Tree integration

`BoundingBallTree` uses `SEBB` to compute node bounding balls. `boundingsphere(tree, node)`
returns a leaf's stored ball unchanged and, for internal nodes, the smallest ball enclosing the
immediate child balls: exact when the solver certifies it, and otherwise the strictly enclosing
approximation described above. `updateradii!` updates nodes bottom-up (children before parents)
so each parent encloses its already-updated children. The two-ball helper
`boundingsphereofspheres` delegates to the same solver.

Only the *containment* half of that contract is load-bearing for the tree. A parent ball has to
enclose its children, and both paths guarantee that by construction; minimality only decides how
tight the result is. So an approximate node ball is slightly larger, which makes slightly more
pairs classify as near — a marginally more expensive near field, never an incorrect partition.
The near predicate's monotonicity argument (see `isnearradius`) likewise rests on containment
alone, not on the radii being minimal.

## Reference

Kaspar Fischer and Bernd Gärtner, *The Smallest Enclosing Ball of Balls: Combinatorial
Structure and Algorithms*, 2004. The relevant results are Lemmas 2.1 (existence/uniqueness),
2.2 (internal tangency + convex-hull center), 2.4 (support of size ``\le d+1``), 2.5 (affine
independence of a basis), and 3.1 (basis ball via a linear system and a quadratic).
