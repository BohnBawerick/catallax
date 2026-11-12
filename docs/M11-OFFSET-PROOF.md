# M11 offset-distance proof route

## Scope and verdict

This is an opt-in experiment for the new M11 proof route. It does not port
blend units or change production fillets. M11 still requires a certified
fillet on an imported spline part and on a loft, with a whole-face rolling-ball
deviation of 100 nm target and 1 micrometre hard cap. The neighbouring-feature
cap and the 1 nm decision tolerance also stay unchanged.

The owner chose the Yang and Jia 2026 engine with an exact per-cell proof
layer. For fillets, apply that engine to the two inward offset surfaces.
The six-variable search in [M27](M27-SPLINE-FILLET-SCOPE.md) remains a recorded
baseline, not the proposed production route. The experiment below tests
whether removing its two normal-length coordinates is sufficient.

It is not sufficient in this experiment. All three edges refuse at depth 24,
inside 20M work, with no certified regular graphs. The 2026-plus-proof route
still needs boundary-aware graph certificates, exact clipping and joins.
No spine, rails or blend face is published, so no whole-face deviation is
certified to either 100 nm or 1 micrometre. This result does not reverse the
owner's engine decision or establish a geometric impossibility.

The spine of all three edges has since been proved by the M10 engine with
exact offset inputs, and a rolling-ball fillet face on each spine; the next
two sections give the results.

## The spine by the M10 engine

Since 12 November the M10 pair search takes exact offsets as input
([offset inputs](M10-SPLINE-SPLINE.md#offset-inputs)). It keeps each
support's parameter rectangle, carries `L = |N|` as an unknown with
`L^2 = N . N`, and certifies a cell only where the offset is proved regular.
Posed as the inward offsets at `r = 1/10` of their exact supports, all three
edges below are proved, each as one open branch of three certified arcs.

| Edge | Answer | Work / 20,000,000 | Search |
| --- | --- | ---: | --- |
| M27 loft edge | 1 open branch, 3 arcs | 2,332,939 | 13 cells, depth 3 |
| Next loft edge | 1 open branch, 3 arcs | 2,666,401 | 10 cells, depth 3 |
| Biquadratic edge | 1 open branch, 3 arcs | 6,793,048 | 13 cells, depth 3 |

Each branch starts on the second support's border `t = 0` and ends on the
first support's border `v = 1`, at exactly isolated roots. On the M27 edge
the ends lie near `(u, v, s, t) = (0.94997, 0.02008, 0.05470, 0)` and
`(0.94181, 1, 0.03853, 0.97520)`. Arcs join at exactly isolated border
points, and a branch ends only on a support border. So the boundary clipping
and the joins that the experiment below left open are proved for the spine. An
independent binary64 Newton solve reproduces every end to six digits. The
records, including every enclosure and the work, are identical at 1, 2 and 4
workers.

```sh
cargo test --release -p catallax-shapes --test offset_pair -- --include-ignored --nocapture
```

## The fillet face on the proven spine

`SplineFillet::fit` builds the rolling-ball face on that spine. A
`SplineFilletDefinition` names two single-span rational patches, the ball's
radius, the side of each support the centre lies on, and a range of one
parameter `tau` of the first support. Each section of the face is the
circular arc round the centre from one contact to the other.

The represented face is a rational B-spline of degree `(D, 2)`:

- Each rail is its support composed with cubic Hermite curves of the contact
  parameters, so it lies exactly on its support. The carrier's boundary rows
  are those compositions, built in exact arithmetic.
- The section's middle row and weight are cubic Hermite fits of the exact
  arc's middle control point and its weight `cos(phi/2) sqrt(W_A W_B)`.

Each span of `tau` proves the following.

- **Spine.** A parametric Krawczyk certificate over the span proves one smooth
  spine point for each `tau`. Interval Taylor jets of that implicit solution
  bound its derivatives to fourth order. Before any span, the M10 search
  proves that the offsets meet in exactly one branch, covering the range.
- **Deviation.** The two-sided position bound follows from the Hermite
  remainders, the rounded node data and the weight errors. Carrier and exact
  face share their parameters, so the bound holds from each surface to the
  other. `SplineFillet::remeasure` recomputes it from the retained node data.
- **Rails on their faces.** Each exact contact's parameters are proved
  inside its support's domain. So is each returned rail, the rounded Hermite
  cubic of each support parameter over each span, which can overshoot
  between nodes that are inside. A cubic lies in the hull of its Bernstein
  coefficients, so a piece whose coefficients are inside is inside, and an
  end coefficient outside is a point of the rail outside. Other pieces are
  halved exactly, up to twelve times.
- **Regular carrier.** The carrier's normal is proved nonzero, from an
  expansion about each box's centre in both parameters.
- **Clearance.** Each contact is proved its support's closest point to the
  centre. The squared distance is convex on a box round the contact, and
  every other cell of the support's domain lies farther than the radius. So
  the ball touches each support face once and enters neither.

The fit starts from eight uniform spans and halves each span whose bound
misses the 100 nm target, at most five times. A fit that meets only the
1 micrometre cap is accepted and reported by its bound. A larger bound
refuses with `fillet.spline-deviation` and the bound reached, in whole
nanometres. Every other failed proof refuses with its own code:
`fillet.spline-spine`, `fillet.spline-rail-domain`,
`fillet.spline-carrier-irregular` or `fillet.spline-clearance`.

| Edge | Range of `v` | Spans | Bound | Work / 20,000,000 |
| --- | --- | ---: | ---: | ---: |
| M27 loft edge | `[1/32, 1]` | 13 | 78 nm | 5,252,648 |
| Next loft edge | `[3/64, 1]` | 16 | 5 nm | 6,970,934 |
| Biquadratic edge | `[1/32, 1]` | 22 | 80 nm | 16,433,004 |

Bounds are rounded up to whole nanometres, and the work includes the spine.
Each range starts just past the spine's end on the second support's border,
near `v = 0.0201`, `0.0407` and `0.0199`. Records, carriers and rails are
identical at 1, 2 and 4 workers and the default. On the M27 face an
independent binary64 solve of the centre puts 297 carrier points on the ball,
and both rails on their contacts, within the proved bound.

**Fitter verdict.** The fit uses Catallax's own fitting method, and it covers
these blends. The existing rolling-ball fitters (`RollingEnvelope`,
`RollingSaddle`, `RollingSphere` and the pipe, torus and mitre faces) each
take one analytic support pair with a closed-form spine, so none accepts a
spline support. Their method carries over unchanged: exact boundary rows on
the supports, a Hermite-fitted middle row and weight, and a bound recomputed
from interval Taylor jets of the exact definition. With the spine from the
M10 engine and jets from a parametric Krawczyk solve, it meets the 100 nm
target on all three edges. No general or external fitter was needed.

## The face between the loft's caps

On the loft each edge runs from the cap `z = 0` to the cap `z = 2`, and
the caps cut the fillet face along curves that are not parameter lines.
On the M27 edge the first rail meets `z = 0` at `tau = 0`, but the second
only at `tau` near 0.0201. Near `tau = 0` a small lobe of the face also
rises above the bottom cap. So the face is fitted past both caps and cut
by them.

- **Extended supports.** `RationalPatch::extended` writes a one-span
  patch's exact polynomials as Bernstein controls on a larger rectangle,
  and proves the new patch like any other. The loft's sides are extended
  by a quarter of the edge each way, to `v` in `[-1/4, 5/4]`, and the face
  is fitted an eighth past each cap, over `tau` in `[-1/8, 9/8]`.
- **Rail borders.** A `RailBorder` names a line where a rail must meet its
  support's border, here `v = 0` and `v = 1` on each side. The fit puts a
  node there. For the first rail the node is that value of `tau` itself.
  For the second it is a binary64 node found along the spine, where the
  rail's parameter is set exactly to the border's value; the two-sided
  bound measures that value against the exact one. So each rail meets each
  cap at a rational point.
- **Cuts.** `SplineFillet::plane_cut` proves where a plane cuts the face.
  On each span, the plane's equation applied to the face is a tensor
  Bernstein polynomial whose coefficients come from the exact rows. Exact
  coefficient signs prove that the face lies on one side of the plane
  before the cut and the other side after it. In between, the face is
  strictly monotone along each curve of constant second parameter. So
  each such curve crosses the plane once, and the cut is a graph
  `tau(theta)`. Its ends must be exact zeros on both rails at fit nodes,
  or the cut refuses with `fillet.spline-cut`.

| Edge | Spans | Bound | Work / 20,000,000 | Bottom cut ends | Top cut ends |
| --- | ---: | ---: | ---: | --- | --- |
| M27 loft edge | 18 | 27 nm | 8,199,428 | `tau` 0 and 0.020075 | `tau` 1 and 1.024669 |
| Next loft edge | 18 | 28 nm | 9,685,010 | `tau` 0 and 0.040749 | `tau` 1 and 1.036694 |

The work includes the spine on the extended supports. Each cut costs about
4,800 units more. Records, carriers, rails and border nodes are identical at
1, 2 and 4 workers and the default. On the M27 face the independent binary64
check covers the whole fitted range.

**The face as a body surface.** `SplineSupport::SplineFillet` holds a fitted
face, like the other fitted rolling-ball faces. Every copy shares what the fit
proved, and a rigid move changes only its frame. A spatial box's chart
preimages come from subdividing the carrier's enclosed Bezier patches. A point
is decided exactly on the face only on a rail, by the rail's exact preimage,
and there its normal is its support's, turned toward the ball. Any other point
is proved off the carrier or refused, as for the saddle and sphere faces. The
carrier is quadratic across, so the line solver of those faces applies: its
span equations come from the exact rows, and each rail's contacts with a line
come from its support's own line hits.

```sh
cargo test --release -p catallax-shapes --test spline_fillet -- --include-ignored --nocapture
```

## The fillet on a loft body

A `fillet` step that selects one edge between two faces on one-span patches,
a patch or a one-chart atlas, rounds it with this face. The edge must run
along a whole border row of both patches, and each of its end vertices must
meet one more face, a plane with rational coefficients, as a loft's caps do.
Each side face must be its whole patch, its four complete border rows and
nothing else, checked from its edges by `RationalPatchFace::of_chart_face`.
So a rail proved inside the patch's rectangle is inside the face. A face an
earlier round trimmed refuses with `fillet.unsupported-spline`.

- **Definition.** The ball rolls inside the material, so each side is the
  opposite of its face's outward normal. Binary64 directions at the edge's
  middle choose the construction: a concave edge refuses with
  `fillet.unsupported-spline`, and the supports are ordered so that the
  carrier's own normal points to the ball. `SplineFillet::faces_ball` then
  proves that order exactly where the first rail starts, from the carrier's
  exact controls and the support's exact normal. Otherwise the step refuses
  with the same code. Both patches are extended by a quarter of the edge each
  way along it, and the face is fitted an eighth past each end, with a rail
  border at each of the edge's two rows on each patch. These are the margins
  of the previous section.
- **Rails on the faces.** Between its two border nodes each rail is proved
  inside its face's own rectangle, from the exact Bernstein coefficients of
  its cubics, or the step refuses with `fillet.spline-rail-domain`.
- **Ends.** Each end plane's cut must end exactly at the two rails' border
  nodes there. The rail points at those nodes are exact controls of the
  carrier. Each must lie strictly inside its side face's straight border
  edge in that plane, by exact collinearity and betweenness.
- **Topology.** Each side face's edge becomes its rail, an exact rational
  B-spline cut from the carrier's boundary row, running the way the edge
  ran. The side faces' border edges in each plane end at the rail points.
  Each plane gains the cut where the edge's vertex was, as an intersection
  edge of the plane and the fillet surface, witnessed at a rail point. The
  blend face lies on the fillet surface, reversed, since the surface's
  normal points to the ball. It is named by the edge's origin pair, so a
  radius edit keeps every face and edge name. The shared edge-blend
  replacement builds the body and checks its topology.

On the M27 loft, the edge between the first two sides is rounded at
`r = 1/10` in 10,349,665 of 20M units, counting the loft, the side faces'
checks, the fit and both cuts. The body has 7 faces and 15 edges. The loft's charts run `u` over
`[0, 1/4]` and `[1/4, 1/2]`, a quarter of the profile each, where the shapes
tests use `[0, 1]`. So the face is fitted on other charts of the same
surfaces: 18 spans and a bound of 59 nm. A radius edit to `3/32` keeps every
name, with 18 spans and 56 nm in 11,075,620 units. Bodies, names and work
are identical at 1, 2 and 4 workers and the default, and so is the refusal
of a second round on a face the first one trimmed.

Not yet done: a point query on this body refuses at a cap's cut edge with
`intersect.parametric-curve`, so Booleans with it wait; STEP into
OpenCASCADE and an imported part remain.

```sh
cargo test --release -p catallax-features --test m11_spline_fillet -- --nocapture
```

The sections below record the earlier experiment unchanged.

## Exact definition and the engine adapter

For support patches A and B with consistently outward normals, retain the
exact definitions

```text
NA = Au cross Av                     NB = Bs cross Bt
OA(u,v) = A(u,v) - r NA / sqrt(NA.NA)
OB(s,t) = B(s,t) - r NB / sqrt(NB.NB)
F(u,v,s,t) = OA(u,v) - OB(s,t)
spine = the selected inward component of F = 0
```

The positive square root selects the normal length. A proved positive lower
bound on each squared length makes this four-parameter definition equivalent
to M27's five equations with two extra positive length variables. It neither
squares away the sign nor fits a spline to the offset. Offset surfaces are
generally not rational spline surfaces, so an ordinary rational-patch input
adapter cannot represent them exactly.

The M10 engine therefore needs an offset evaluator that returns rational
enclosures of position and derivatives, with its original exact support
definition attached. For a support parameter i, the derivative is

```text
Di OA = Di A - r [Di NA / L - NA (NA.Di NA) / L^3]
L = sqrt(NA.NA) > 0
```

The experiment encloses each positive square root by 32 exact rational
bisections. These intervals enclose the exact definition; their width is not
a modelling tolerance. A denominator that cannot be separated from zero
returns `experiment.m11-normal-range-unproved`. That result means the bound
failed, not that the support has a proved singularity. A production adapter
must also prove offset regularity and the selected side. A regular support
alone does not exclude a focal singularity of its offset.

## Distance-field candidates and complete proof

Yang and Jia's [2026 paper](https://doi.org/10.1145/3807948), equations 2 and 3,
uses the oriented distance to the unique closest point on the second surface.
For offsets, replace those surfaces by OA and OB. The field components are
the target offset unit normal dotted with the two source offset derivatives.
A positive common normalization factor does not change the field winding.

The experimental candidate stage uses a 3 by 3 target grid. It chooses the
nearest grid point by the midpoint of its exact squared-distance enclosure,
with fixed index order for ties. It evaluates the two field components there,
classifies their quadrants and counts the cyclic quadrant changes. An edge
joining opposite or unresolved quadrants causes four-way source subdivision,
up to two levels. Every new source sample repeats the nearest-grid search.
Ambiguous nearest-grid choices are counted. This implements the grid proposal
part of the method, not a proof that the selected point is the continuous
closest point.

Neither a zero candidate winding nor exhausted candidate subdivision proves
an empty cell. Section 6 of the paper explicitly limits grid coverage and
discusses multiple critical points within one cell. Catallax's exact layer
must cover the full pair of offset domains independently of those candidates.
Before a winding can become a certificate, it needs continuous closest-point
uniqueness, offset regularity, a nonzero boundary field and a proof that no
unclassified critical components remain. A sampled nearest point supplies
none of those premises. Non-isolated critical curves need separate treatment.

The experiment searches all of `[0,1]^4`, without a corridor around the known
edge. Rational interval evaluation of the exact support polynomials and
normalized offsets can exclude a box only when some component of F has a
strict sign. Otherwise it tries the same regular-graph certificate as M27,
now with a three-by-three solved Jacobian. For each free coordinate in a fixed
order, a rational midpoint inverse P gives

```text
T(Y,z) = m - P F(m,z) + (I - P Jy(Y,z)) (Y-m)
```

Strict inclusion in Y and a strict row-sum contraction bound below one prove
one regular equality graph across the whole free interval. Failure causes
midpoint subdivision of the longest parameter interval, with axis order
breaking ties. All arithmetic reservations use the existing numeric work
prices before execution. Coefficient scans, interval operations, candidate
evaluation, inversions, graph tests and subdivision are charged to one shared
20,000,000-unit request. Immutable symbolic fixture construction and the
independent control checks are outside that count, as in M27.

The caps remain depth 24 and 16,384 proof cells. A remaining box at either
cap refuses with `experiment.m11-offset-depth-limit` or
`experiment.m11-offset-cell-limit`. Numeric exhaustion retains the existing
`numbers.operation-budget`. No timeout is a geometric outcome.

## Boundary witnesses, rails and the fillet face

A complete cover of equality graphs is only the first certificate. It must
retain exact boundary roots, shared endpoint witnesses and branch joins.
Graphs which touch or leave a support domain need proved clipping. Strict
interior interval inclusion cannot by itself certify an arc through a solved
coordinate boundary. Until these exist, a nonempty completed cover refuses
with `experiment.m11-offset-joins-unproved`.

On a certified, selected center branch C, construct the rails directly:

```text
RA = A(u,v)                         RB = B(s,t)
RA - C = r NA / sqrt(NA.NA)         RB - C = r NB / sqrt(NB.NB)
```

These identities prove radius and support tangency. Do not send the rails
back through a crossing-intersection classifier, since their tangency is a
construction fact. Retain their support parameters and common center witness.
A normal contact alone does not prove an empty rolling ball. The selected
rail must also be a closest point on its support, with the material side and
clearance from other faces proved.

M11 must then select a consistent circular arc and prove that the fitted
blend covers that rolling-ball definition in both directions. The face,
rails, ends, material side and absence of folds all need certificates.
Report the reached whole-face deviation if it exceeds the 100 nm target or
1 micrometre cap. Offset equality, local radius residuals and sampled points
cannot substitute for that deviation. Shared trims, valid solid topology,
stable names under a radius edit and independent STEP import remain later
M11 gates. Production singular-offset, ambiguous-projection, collapsed-radius
and incomplete-join refusals must remain distinguishable.

## Public cases

All cases use radius `r=1/10` mm and unit rational spline weights.

1. `m27-loft-edge` is exactly M27's first two bilinear side patches and full
   parameter domains. No affine surrogate replaces that loft edge.
2. `loft-next-edge` uses the same loft's second and third side patches. The
   third patch is `C(s,t)=(1-2s-t, 1+t(1-2s), 2t)`. Its shared edge with B is
   `B(1,h)=C(0,h)=(1-h,1+h,2h)`.
3. `biquadratic-edge` adds `u(1-u)v(1-v)/5` to M27's A.z and subtracts
   `s(1-s)t(1-t)/5` from B.z. Both curved spline patches retain the original
   exact shared edge. This synthetic higher-degree control is not an imported
   part or evidence of arbitrary imported-face acceptance.

## Measured results

Each row returns `experiment.m11-offset-depth-limit` at depth 24. None exhausts
the allowance; all certify zero regular graph cells. Candidate winding is
zero on the initial source cell in every case, with four source samples,
no unresolved candidate edge and no ambiguous nearest-grid choice. The exact
layer still searches the full domain. A zero grid index never skips that work.

| Edge | Work / 20,000,000 | Candidate work, included | Visited / empty / pending cells | Last box offset-gap upper bound |
| --- | --- | --- | --- | --- |
| M27 loft edge | 5,431,053 | 47,321 | 195 / 91 / 13 | 148,170 nm |
| Next loft edge | 4,518,992 | 48,939 | 159 / 71 / 17 | 147,813 nm |
| Biquadratic edge | 10,168,139 | 78,593 | 209 / 98 / 13 | 151,770 nm |

The last column is a rational interval upper bound on
`abs(Fx)+abs(Fy)+abs(Fz)` over the first unresolved box, rounded upward to a
whole nanometre. It bounds offset separation for every parameter pair in that
box. It is not a fitted-face deviation or a lower bound on attainable error.
The reached whole-face deviation is **unproved** for all three rows; neither
the 100 nm target nor the 1,000 nm hard cap has passed. Reporting these local
gap bounds as fillet accuracy would be incorrect.

The first and third rows stop in the box

```text
u in [59/64,15/16]     v in [7/64,1/8]
s in [3/64,1/16]      t in [3/32,7/64]
```

The second row stops in

```text
u in [59/64,15/16]     v in [1/32,3/64]
s in [3/64,1/16]      t in [0,1/64]
```

All parameter widths are 1/64. Positive normal-length enclosures remove two
search dimensions, but this certificate still cannot finish the branch cover.
The fresh M27 baseline also reproduces its published 2,910,656 work units,
102 visits, 44 empty cells, zero graphs and depth refusal at 1/2/4 workers.
These are different proof implementations and work totals, not a speed
comparison or a measurement isolating the effect of dimensional reduction.

The M11 run passes complete-record comparison at 1, 2 and 4 workers. Small
controls pass as well. Published work, cell counts, candidate evidence, refusal
code and diagnostic bound ceilings are checked by the opt-in test.

The immediate next proof work is to handle solved-coordinate boundaries and
attach exact root and join witnesses to accepted cells. After that, test
continuous closest-point uniqueness and offset regularity before promoting
winding evidence to a certificate. The blend port and imported-part gate
remain open; increasing the allowance or changing a tolerance would not supply
these missing proofs.

## Reproduction and checks

```sh
cargo test --release -p catallax-shapes --test m11_offset_spine -j2 -- --include-ignored --nocapture
cargo clippy --release -p catallax-shapes --test m11_offset_spine -j2 -- -D warnings
cargo fmt --all -- --check
```

The scope run is ignored by ordinary tests and the scoreboard. Small controls
check square-root enclosure, signed squaring, cyclic winding including the
fourth-to-first quadrant transition, the exact M27 boundary offset and its
derivatives, a regular equality graph, rank deficiency and budget exhaustion.
Independent requests compare complete ordered records at 1, 2 and 4 workers,
including work, candidate evidence, terminal refusal and unresolved box.
This schedules cases across workers; it does not qualify parallel work within
one offset search or claim a speed pass.

Local execution used `rustc -O --test` against the Cargo-built release numeric
library. Both M11 tests and the unchanged M27 baseline passed. The packaged
release Clippy target passed with warnings denied, and workspace formatting
passed. The full local Cargo test build stopped at its 900-second build
timeout during shapes code generation; that incomplete build is not a
geometric refusal or a packaged test pass. The 477-case replay was not run
locally for this test-and-documentation change.

At the coordination checkpoint, the M10 spline-spline lane had no pushed head
to read. This experiment uses the landed R13 equations and the supplied 2026
paper. Integration with M10's offset evaluator and shared proof objects
remains open; this change edits none of the M10 lane's files.

No production code or replay input changes. The 477 saved replay answers are
untouched. Raw evidence retained privately.
