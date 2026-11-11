# M27 spline-fillet proof scope

## Verdict

Change the proof route before committing M11's arbitrary spline-face scope.
The direct six-variable extension of M10 tested here refuses at proof depth
24 after 2,910,656 of 20,000,000 application work units. It does not exhaust
the work allowance. This is evidence against that particular
search, not a proof that every exact offset search is unaffordable. M11's
roadmap text stays unchanged for the owner's later decision. No blend unit is
ported by this scope.

## What M10 already proves

[The plane slice](M10-PLANE-SECTIONS.md) and
[the cylinder and sphere slice](M10-QUADRIC-SECTIONS.md) substitute an analytic
surface equation into one rational patch. There is one polynomial in two
parameters. Exact Bernstein signs exclude empty cells; a nonzero partial
and two simple boundary crossings certify a regular arc. Exact boundary
roots and endpoint joins turn the full cell cover into components. The
search starts without a tracer and finds the small loops a tracer misses.

That last property is R1's requirement. R2 leaves singular contact as a typed
refusal until exact contact isolation can classify it. R13 asks separately
how to cover spline-against-spline intersections in four parameters. None of
those slices supplies a general five-equation, six-variable curve search.
The present experiment extends the cell search in test code only. It does
not claim that the production M10 entry accepts offset surfaces.

## The loft edge and its offset equations

The case is `nonaffine_parallel_sections_use_a_whole_polar_embedding_proof`
in `crates/catallax-app/tests/rational_lofts.rs`, listed in
[LOFT-CASES.md](LOFT-CASES.md). The bottom polygon at height zero has corners
`(-1,-1), (1,-1), (1,1), (-1,1)`. The top polygon at height two has corners
`(0,-1), (1,0), (0,2), (-2,0)`. The ruled correspondence gives these first
two side patches, with every parameter in `[0,1]`:

```text
A(u,v) = (-1 + 2u + v - uv, -1 + uv, 2v)
B(s,t) = (1 - st, -1 + 2s + t, 2t)
```

Both are nonplanar bilinear spline faces, with unit rational weights. Their
common edge is `A(1,h) = B(0,h) = (1,-1+h,2h)`. This tests an edge between
two spline faces, rather than a spline against a planar cap. The experiment
reconstructs these charts exactly from the existing case; it does not build
the application body, display, trim topology or STEP.

Their outward, unnormalised normals are the exact cross products:

```text
NA = Au x Av = (2v, 2v-4, 2u-v)
NB = Bs x Bt = (4, 2t, 2s-t)
```

Choose an inward rolling ball of radius `r = 1/10` millimetre. Its center
must satisfy `C = A - r NA/a = B - r NB/b`. Keep the square roots as two
extra unknowns, `a = sqrt(NA.NA)` and `b = sqrt(NB.NB)`. The five polynomial
equations in the ordered variables `(u,v,s,t,a,b)` are:

```text
10 (Ai - Bi) a b - NAi b + NBi a = 0, i = x,y,z
a^2 - NA.NA = 0
b^2 - NB.NB = 0
```

The complete finite search box is `[0,1]^4 x [2,5] x [4,5]`.
Throughout the two patches, `NA.NA >= 4`, `NB.NB >= 16`, and each is at
most 24, by the displayed component bounds. Thus this box contains every
positive normal length. Neither zero nor the negative square-root sheets
can enter. Multiplication by `ab` loses no solution and adds none in this
box. The entire pair of side domains is searched, including possible remote
components, with no guessed corridor around the edge.

The opt-in test checks the charts against bilinear interpolation of their
corners, the normal polynomial identities against the two cross products,
and the cleared equations against the original offset differences.

## Search and proof limits

The experiment keeps M10's exact tensor Bernstein exclusion and midpoint
de Casteljau subdivision, depth-first with the left child first. It splits
the longest absolute coordinate interval, breaking ties by variable order.
The unchanged limits are depth 24, 16,384 proof cells and 20,000,000 operation
units. It never substitutes a wall-clock timeout for deterministic work.

M10's scalar derivative test cannot certify this system. For each possible
free coordinate, in a fixed order, the experiment tries a rational
preconditioner for the other five variables. It bounds their Jacobian with
Bernstein derivative coefficients. For fixed free coordinate `z`, let
`F(y,z)` be the five equations, `m` the solved-box midpoint and `P` an
invertible rational midpoint-Jacobian preconditioner. The interval image is:

```text
T(Y,z) = m - P F(m,z) + (I - P Jy(Y,z)) (Y-m)
```

Strict inclusion of this image in `Y`, together with a strict row-sum
contraction bound below one, proves exactly one regular graph across the
whole free interval. Otherwise the cell stays unresolved and splits. The
residual range uses all free-coordinate Bernstein coefficients after
restricting the other coordinates to their midpoints. It is not a sampled
residual or a seed-dependent tracing test.

Every proof addition, multiplication and division reserves the existing
numeric exact-arithmetic charge before executing. Explicit checkpoints
charge coefficient scans and cell processing. The measured request includes
Bernstein conversion, exclusion, derivative bounds, rational matrix
inversion, graph tests and splitting. The small symbolic fixture definition
and its independent self-checks are not an application loft build and are
not included in the search count. No tolerance or work-price table changes.

A completed regular-cell cover would still need exact boundary witnesses
and joins, as M10 does. The experiment refuses that unimplemented last step
instead of publishing disconnected graph boxes as a spine. A resource
refusal leaves all remaining cells unresolved; neither empty geometry nor a
partial fillet is returned.

## Measured result

| Measure | Result |
| --- | --- |
| Finished the complete search | No |
| Final typed refusal | `experiment.m27-depth-limit` |
| Cause | Neither Bernstein exclusion nor a regular graph was proved at depth 24 |
| Accepted operation work | 2,910,656 / 20,000,000 units, 14.55328% |
| Initial Bernstein conversion | 24,550 units, included in the total |
| Work allowance exhausted | No, 17,089,344 units remained |
| Visited cells | 102 |
| Proved empty cells | 44 |
| Certified regular graph cells | 0 |
| Unresolved cells at refusal | 14, including the cell being processed |
| Deepest binary split | 24 |
| Published spine, rails or blend | None |

The first unresolved box at that depth, in the declared variable order, is:

```text
u in [7/8, 15/16]       v in [7/16, 1/2]
s in [0, 1/16]         t in [3/8, 7/16]
a in [53/16, 109/32]   b in [4, 33/8]
```

All seven requests, one serial and batches of two and four concurrent
requests, return the same refusal, box, work and cell counts. There is one
terminal refusal per request. No resource, tangent or join refusal occurs
in this run. The experiment also has `experiment.m27-cell-limit` for the
unchanged cell cap, `experiment.m27-joins-unproved` if a nonempty regular
cover finishes without exact joins, and the numeric layer's existing work
refusals. Those paths are not observed results for this edge.

This is a depth/coverage failure, not evidence that the spine needs more
than 20M units. The remaining work to finish it is unknown. At depth 24 the
six coordinate intervals are still roughly 1/16 to 1/8 wide; distributing
splits among six variables is different from the plane slice's two. The
regular-graph test is sufficient and more restrictive than M10's scalar
two-boundary-crossing test. Its failure is not a geometric singularity
certificate and does not establish that a rolling-ball fillet is impossible.

The test compares complete result records from independent requests at 1,
2 and 4 workers. This measures answer determinism under concurrent requests;
it does not claim a parallel speedup or a shared multiworker work allowance.
Raw evidence retained privately.

## What a fillet proof must add

The two contact rails are built from the center branch's support parameters:
`RA = A(u,v)` and `RB = B(s,t)`. Their radius vectors are `r NA/a` and
`r NB/b`, so radius and support tangency follow from the retained definition.
Running a crossing intersection solver on the rails would ask it to classify
an everywhere-tangent contact and would discard that construction fact.

A publishable fillet still needs all of the following:

- A complete center-branch cover with exact shared endpoint witnesses,
  trim-domain checks and a selected inward side. Singular offset crossings,
  support singularities and radius collapse need distinct typed refusals.
- A consistent circular arc between the two radius vectors at every center,
  with the chosen material side retained. A local normal-length proof alone
  does not select that arc or exclude a folded or self-intersecting blend.
- A whole-face, two-sided distance certificate between any fitted blend and
  this rolling-ball definition, including end caps and both rails. The target
  stays 100 nm, the hard cap 1 micrometre, and the local feature-size rule
  stays one tenth of the shortest neighbouring feature. A failed check must
  report the bound reached, rather than accepting a fit from point samples.
- Shared rail objects for neighbouring trims, checked solid topology, stable
  names under a radius edit and a STEP round trip through OpenCASCADE. These
  belong to M11 and are not supplied by a centerline experiment.

## Consequence for M11

The next proof experiment should add boundary-aware graph certification
and exact endpoint joins, and test dimension reduction or tighter global
normal and correspondence bounds before local graph certification. Candidate curves may propose local
boxes, but the exact global search must still exclude every other component.
An exact elimination or a problem-specific parameter reduction must retain
positive normal lengths and prove equivalence to the original offset system.
An approximate offset fit alone does not meet R1 or the rolling-ball gate.

Keep the arbitrary spline blend port behind that proof experiment. It must
finish the same edge, include exact joins, publish work inside the unchanged
20M allowance and handle a singular-contact control with a typed result.
Then test a genuinely higher-degree or imported pair before treating this
bilinear case as evidence for arbitrary spline faces. The present refusal
supports changing the proof route, not raising the allowance, narrowing a
local tolerance, or silently dropping small components. It does not settle
R13's separate spline-against-spline search decision.

## M11 offset-distance follow-up

The next experiment is now [the M11 offset-distance route](M11-OFFSET-PROOF.md).
It applies the owner's 2026-engine decision to the intersection of the exact
inward offsets, with positive normal lengths enclosed rather than searched
as independent coordinates. This changes the search from six parameters to
four. The 2026 distance-field winding and subdivision propose starts; exact
per-cell exclusion and graph checks retain authority over complete coverage.
Contact rails follow from the center branch's retained support parameters.

The M27 result above remains the six-variable baseline. The new experiment
uses the same loft edge, its next spline-spline edge and a biquadratic control,
with the unchanged 20M allowance. Its note records the verdict, refusal,
uncertainty and work for each case, including 1/2/4-worker comparison. It is
not a production fillet port or a completed M11 gate.

All three new cases refuse at depth 24 with no regular graphs. The original
loft edge uses 5,431,053 of 20M units; its last unresolved box has parameter
widths of 1/64. No whole-face deviation is certified to the 100 nm target or
1 micrometre cap. Removing the two length coordinates alone does not finish
the proof. Boundary-aware graphs, clipping and exact joins remain necessary.

## Reproduction

```sh
cargo test --release -p catallax-shapes --test m27_spline_fillet -- --ignored --nocapture
```

The test is ignored by default. Ordinary tests and the scoreboard do not run
it. The original scope added the M27 ledger line in `PLAN.md` and kept its
search entirely in a new integration test. The M11 follow-up adds its own
opt-in test and proof note. No production geometry behavior, saved replay
answer, allowance or tolerance changes.
