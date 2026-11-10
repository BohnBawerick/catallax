# M10 plane sections

The first slice adds exact regular branches for a single `RationalPatch` cut
by a plane. Existing univariate and affine reductions still return their exact
rational splines. A general regular cut returns `Curve::PatchPlane`, with the
patch, plane and ordered parameter-cell cover retained as its definition.
The fitted controls never define the meeting geometry.

## Complete search

Each knot span starts with the exact Bernstein coefficients of the plane
equation after homogeneous substitution. Strictly single-sign coefficients
prove a cell empty. A strictly single-sign partial derivative and exactly two
simple boundary crossings prove one regular arc. Every remaining cell splits
along its longer parameter direction, under the existing depth, cell and work
allowances. Exact edge roots join the arcs across cells and knot lines. A
vertex has degree two, or degree one on the exact patch border. Traversing
that graph gives every open component and closed loop.

This search runs without a candidate. It finds the hill's loop one millionth
of a millimetre below the peak even when the pinned tracer returns no curve.
At the exact hilltop it refuses with `intersect.patch-plane-tangent` and names
`(5/2, 5/2, 225/256)`, at parameters `(3/2, 3/2)`. Other unresolved contacts
and exhausted allowances remain typed refusals.

Membership uses exact patch inversion and the component's cells. Open order
and closed cyclic order follow the cell traversal. Tangents use the implicit
parameter equation, retaining the arriving and leaving knot sides. Pair
orientation follows the plane-normal cross patch-normal direction and
reverses when the input surfaces swap.

## Candidate certificate

`PatchPlaneBranch::certify_candidate` checks an open single-span graph branch.
The representation is the supplied cubic parameter curve composed with the
retained rational patch. Each piece must be strictly monotone in the free
parameter, join its neighbour exactly and reach both exact patch borders.
The selected branch endpoints must lie in the tube's end caps.

For each piece, exact Bernstein signs prove opposite cut-equation signs on
the two transverse walls. A positive transverse derivative lower bound over
the full tube rectangle proves uniqueness and bounds the transverse inverse.
The exact branch follows the tube between its end caps. Monotone free-parameter
coverage gives correspondence in both directions, including endpoints.
A rational upper bound on surface speed converts the tube half-width `rho`
to a whole-component Euclidean distance bound `U` in millimetres.

The immutable `PlaneBranchApproximation` retains the exact definition,
unchanged candidate controls, `rho`, speed bound, positive slope lower bound,
`U`, recorded cap, work and `GeometricProof::Approximation`. Its proof names
resolve within that record to the retained definition and checked tube theorem.
It accepts the existing 100 nm target or 1 micrometre hard cap. Failure proves
no admissible tube by this checker, so it refuses rather than asserting an
excessive-distance lower bound.

Exact arithmetic reservations use the existing numeric work-cost model and
caller operation allowance. Tests put patch admission, search and candidate
checking inside the unchanged 20M operation allowance. No tolerance or work
limit changes.

## Verification

```sh
cargo test -p catallax-intersect --test patch_plane -- --nocapture
```

The regression cases include the rational loft's exact endpoints, the tight
and default C++ candidates, a candidate shifted by `1e-3`, an end cap that
stops short, two wavy open branches, the hill loop and tiny loop, exact tangent
refusal, an empty cut above the hilltop and work exhaustion. Fixture provenance
and regeneration are in the [oracle fixture notes](../crates/catallax-intersect/tests/fixtures/patch-plane/README.md).
Raw evidence retained privately.

The measured regression operations, including patch admission, are below.
The allowance remains 20,000,000 units for every row.

| Case | Result | Operation units |
| --- | --- | ---: |
| Tilted loft | 1 open branch, 2 arc cells, depth 1 | 96,824 |
| Tight candidate | U < 91.871 nm, rho = 1/65536 | 560,119 |
| Default candidate | Typed cap refusal | 297,195 |
| Candidate shifted by 1e-3 | Typed cap refusal | 97,862 |
| Candidate stopping short | Typed endpoint refusal | 96,836 |
| Wavy surface | 2 open branches, 6 arc and 4 empty cells | 391,253 |
| Hill at z=0.3 | 1 closed loop, 10 arc and 15 empty cells | 152,036 |
| Hill 1e-6 below its peak | 1 closed loop, 12 arc and 76 empty cells, depth 21 | 174,443 |
| Exact hilltop | Typed tangent refusal naming the touch point | 60,998 |
| Hill 1e-6 above its peak | Proved empty | 63,361 |

The 1-, 2- and 4-worker checks agree on exact results, endpoint enclosures,
refusal details, candidate bound and operation units. Their extra endpoint
enclosures give 560,215 units for the loft with the tight candidate, 395,982
for the wavy patch, 153,252 for the hill and 175,755 for the tiny loop.

This slice does not add cylinder or sphere cuts, third-surface root queries,
atlas joins, spline-against-spline intersections, Boolean consumption, or a
STEP writer for these branches. Multi-span and closed branches retain their
exact definitions; the first candidate representation checker refuses them.
