# M10 spline-with-spline meetings

This stage meets two rational spline patches. `RationalPatch::patch_meeting`
returns every branch of the meeting as an ordered chain of certified arcs.
When it cannot prove the whole meeting, it returns a typed refusal that names
the cell where the proof stopped. The proof route is the owner's decision of
11 November 2026: the 2026 field method of Yang and Jia as the engine, with
an exact proof layer under every cell. [R13](R13-SPLINE-TOPOLOGY.md) records
the decision beside the experiment that preceded it.

## Method

The meeting is the solution set of `A(u,v) = B(s,t)` in four unknowns. The
2026 paper follows the field `(Au . Nb, Av . Nb)` on the first patch's
parameters and counts its winding numbers on a grid. Its Section 6 limits
the guarantee to that grid's resolution. Here subdivision replaces the grid,
and a cell is accepted only by an exact certificate, so no fixed spacing
limits what the search finds.

The search splits the first patch's knot spans into cells. Each cell carries
a cover: the second patch's cells, split until they are no larger than the
cell, that exclusion has not removed. The meeting is that of the closed
patches, borders included, so exclusion must prove that the closed cells
hold no meeting point. It does when the boxes round the two control nets are
disjoint, or their slabs along either cell's normal direction are disjoint.

Boxes that touch exactly prove less. Every meeting point then lies on the
touching side, and where one patch reaches that side on one border of its
cell only, on that border. Each patch's cells own their points half-open: a
cell owns its low borders, and leaves its high borders to the neighbour,
except on the physical border, which it owns. A touching pair is dropped
only when its meeting points must lie on an internal high border. They then
belong to the neighbouring cell across it, and a pair that owns a meeting
point is never dropped. A meeting on a physical border is never dropped this
way.

A first-patch cell against the cover cells of one second-patch knot span is a
unit. A unit is accepted only when one of three certificates holds over all
of it:

| Certificate | Strict sign of | Consequence |
| --- | --- | --- |
| `FieldFirst` | one component of `(Au . Nb, Av . Nb)` | `u` or `v` is strictly monotone along every branch in the unit |
| `FieldSecond` | one component of `(Bs . Na, Bt . Na)` | `s` or `t` is strictly monotone along every branch in the unit |
| `NormalCone` | one coordinate of `Na x Nb` | that spatial coordinate is strictly monotone along every branch |

Each certificate rules out a closed loop, an isolated contact and a critical
point inside the unit. For the two field certificates this is the paper's
winding-number test in cell form: a field with one component of strict sign
has winding number zero round every closed curve in the unit. Tangents and
normals are homogeneous, positive multiples of the true ones, so no square
root or division enters. Their bounds come from exact Bernstein nets of each
cell.

A certified unit must also show either no border points or exactly two, one
where the branch enters and one where it leaves. A border point is a root of
the face system: the three equations `A_k B_w - B_k A_w = 0` with one
parameter fixed on a face of the unit. An exact Krawczyk certificate isolates
each root on its face box. A second certificate then tries a small box round
a binary64 Newton estimate, so later sign tests settle cheaply. Each
certificate is capped at
50,000 units. A capped or undecided face splits the cell; it is never
accepted. Where both patches share a chart, a face can meet the other patch
along one of its corners. When an exact combination of the face equations
vanishes on that corner, with a quotient of one strict sign over the face,
every root lies on the corner and the other two equations decide it.

The tangent of a branch is `Na x Nb`. Its parameter signs follow from the two
fields. A border point where the tangent points into the unit on every face
through the point is a start, and one where it points out on every such face
is an end. Any other border point only touches the unit's corner and belongs
to no arc of that unit.

An unsettled cell splits into thirds along its longer side. Second-patch
cells split at two and four fifths. The two rules never share a split value
and never split at a half, so the cell faces of the two patches meet only on
knot lines. Before a split, the search tests the simplest rationals in the
unit for an exact tangent contact, where `A = B` and the normals are
parallel. A contact refuses with `intersect.spline-pair-tangent` and names the
point. A face whose three equations are exactly dependent has a curve of
roots or none. Unless one equation keeps a strict sign over the face, a
branch may run along that cell or patch border, which this stage does not
certify, and the face refuses with `intersect.spline-pair-boundary-overlap`.
A cell still unsettled after 48 splits refuses with
`intersect.spline-pair-unresolved` and names the cell and its cover. No
near-enough answer is ever accepted, and an empty answer proves that the
closed patches do not meet.

Arcs join at shared border points. Two ends are the same point when one's
small box lies inside the other's uniqueness box on the same face, or when one
Krawczyk certificate on their hull proves a single root. Every join point
needs one start and one end, or a single start or end on a physical patch
border. Anything else refuses with `intersect.spline-pair-join`.

Before the search, one exact rotation turns both patches. A binary64 search
over small integer quaternions picks the rotation that best brings the first
patch's mean normal to `z` and its mean `u` tangent to `x`. The identity stays
unless another rotation is clearly better. Parameters do not change under a
rotation, so neither does the answer, and boxes and slabs then fit a turned
copy as closely as the original.

Binary64 arithmetic only chooses. It chooses the rotation, the small box that
the second root certificate then checks, and which exact slab tests to skip
because they would clearly fail. Every accepted fact is exact. Patches with
irrational coefficients refuse with `intersect.spline-pair-irrational`.

## Threads and determinism

The search runs level by level. Before each level, the request's own thread
refines every cell's cover in a fixed order. It builds shared second-patch
children, tangent and normal nets, and span-pair equations once, so their
charges fall in that order. Each cell's proof then runs on helper threads,
one per available core. The request's thread adopts their recorded charges
in cell order, or proves a cell itself when no helper has finished it. The
answer, the refusal text and the work count are those of one thread.

## Verification

```sh
cargo test --release -p catallax-shapes --test patch_pair -- --nocapture
```

The inputs are R13's 24 pairs: the M25 wavy pair and the hill families
`loop`, `contact`, `small-loop` and `disjoint`, plus the open control. Each
comes in a base, moved, turned and scaled copy, with the controls posed
directly; the search sees only the posed patches. R13 gives the constructions.

The tiny loop is the hill pair at `c = 225/256 - 1/10^10`. Its loop has a
parameter diameter of about `2.4e-5`. The paper's general experiment, at
density 15 and depth 1, and R13's 5 by 5 grid at depth 2 both leave cells
about a tenth of a unit wide, thousands of times wider than the loop. The
search finds it as one closed branch in every pose.

The allowance is the unchanged 20,000,000 operation units for every row.
Work includes setup, the search, root certificates and the join.

| Case | Result | Search | Base | Moved | Turned | Scaled |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| M25 wavy pair | 2 open branches | 63 cells, 17 arcs, depth 5 | 3,520,954 | 3,580,394 | 3,531,451 | 3,635,304 |
| Loop, `c = 3/10` | 1 loop | 351 cells, 88 arcs, depth 8 | 16,603,615 | 16,667,992 | 16,610,835 | 16,632,719 |
| Contact, `c = 225/256` | Tangent refusal | refused after 2 splits | 902,975 | 906,266 | 910,195 | 902,975 |
| Small loop, `c = 225/256 - 1/10^6` | 1 loop | 363 cells, 51 arcs, depth 16 | 16,766,542 | 16,852,045 | 16,773,762 | 16,934,205 |
| Disjoint, `c = 225/256 + 1/10^6` | Proved empty | 147 cells, 0 arcs, depth 12 | 3,456,618 | 3,486,561 | 3,463,838 | 3,492,516 |
| Open control | 1 open branch | 12 cells, 3 arcs, depth 1 | 1,020,750 | 1,026,329 | 1,027,970 | 1,021,632 |
| Tiny loop, `c = 225/256 - 1/10^10` | 1 loop | 372 cells, 26 arcs, depth 24 | 16,404,984 | 16,495,755 | 16,412,204 | 16,608,773 |

Cells are the first patch's cells examined, and arcs are certified arcs.
Every row stays inside the allowance; the largest is the scaled small loop at
16,934,205 units. The contact refusal names the exact hilltop at parameters
`(3/2, 3/2)` on both patches and at the posed point in each copy, such as
`(141/64, 5/2, -377/256)` in the turned one. The wavy pair gives two open
branches, matching the two curves SMLib reports. Copies keep the same cells,
arcs and depth. Their work differs slightly because posed coefficients have
more limbs.

Another test runs the six base copies, the turned loop and the tiny loop at
1, 2 and 4 workers and at the default, which is every core. Branches, arcs,
point enclosures, refusal text and work agree exactly. On the loop, every
arc carries one of the three certificates, consecutive arcs share their join
point's enclosure, and every enclosure allows `u = s` and `v = t`, as the
shared chart requires. The loop with every first-patch weight set to 2, the
same surface with a weight that is not one, gives the same arcs in the same
cells. A one-unit allowance refuses with `numbers.operation-budget`.

Two border controls guard the closed-domain rule. In both,
`A(u,v) = (u,v,0)` on `[0,1]^2`, and the boxes touch exactly at `x = 1`,
which `A` reaches only on its physical border `u = 1`. With
`B(s,t) = (1+s, t, s)`, the patches meet along the whole segment `(1,v,0)`,
on a physical border of each, with transverse normals. That pair refuses
with `intersect.spline-pair-boundary-overlap` at the face `u = 1`, in 9,465
units. With `B(s,t) = (1+s+t², t, s)`, they meet only at `(1,0,0)`, a
corner of both domains, where the face system has a double root. That pair
refuses with `intersect.spline-pair-unresolved` at the corner cell after 48
splits, in 3,139,030 units. Before the closed-domain rule, both returned an
empty answer.

## Offset inputs

`RationalPatch::offset(d)` names the exact offset `P + d N/|N|` of a patch,
with `N = P_u x P_v` and a signed rational distance `d`.
`OffsetPatch::offset_meeting` meets two offsets with the same search, over
both supports' closed parameter domains, and reports the supports'
parameters. A rational patch's search is unchanged: its answers and work
equal those in the table above.

An offset is not a rational surface, because `|N|` is the square root of a
polynomial. The search carries `L = |N|` as one more unknown per offset side,
with the equation `L^2 = N . N`. Over a cell, `L` lies in a positive rational
enclosure: the exact Bernstein range of `N . N`, with square-root bounds that
binary64 proposes and exact squaring checks. Each part of the search then
reads the offset as follows.

- **Exclusion.** Over a cell the offset lies in its support's control hull
  moved by `d N/|N|`. The ranges of `N` and `L` bound that move, so each
  control box grows by that bound. The exact touch rule reads exact controls,
  so it applies to rational patches only.
- **Certificates.** Each tangent of the offset, times `W^2 L^3`, is
  `L p1 + p0`, with `p1` and `p0` polynomials, so cells bound it from exact
  nets. A unit is certified only where the triple product of the offset's two
  tangents with `N` is positive on every cell. The offset is then regular
  there, its normal is a positive multiple of `N`, and the field and
  normal-cone certificates read `N`.
- **Border points.** The offset's homogeneous coordinates are
  `(L P + d W N, W L)`, so a face system gains one unknown and one equation
  per offset side. Its box for `L` encloses `L` over the cells, widened on
  both sides by half its width, so that every root lies strictly inside. One
  Krawczyk step on such a system can cost more than the 50,000-unit cap of a
  plain face certificate. Its cap is sixteen of its own steps instead, and
  the request's allowance still pays every unit.
- **Contacts.** An offset point at rational parameters is not rational, so
  offsets have no exact tangent-contact test. An offset pair with a tangent
  contact refuses at the split limit, as does a branch through a corner that
  both supports' borders share.

Two plane offsets meet on their exact line in one arc and one cell, in
101,770 units. [M11](M11-OFFSET-PROOF.md) records the fillet spines of its
three edges.

```sh
cargo test --release -p catallax-shapes --test offset_pair -- --include-ignored --nocapture
```

## Joined atlas inputs

`RationalAtlas::patch_meeting` feeds every retained chart into the same
Yang and Jia cell-pair engine. Each span carries its chart identity and
physical-border flags. An exactly identified periodic high row leaves its
points to the low row, while a physical border still owns its points.

The join stage canonicalizes fixed periodic row coordinates only after the
atlas has proved exact position and G1 gluing. It then checks the same root
uniqueness boxes or certifies their hull. The returned root parameters keep
their original charts; only the graph's join keys use the identification.
`PairArc` records both chart indices and each chart-local span index.

`RationalAtlas::partition` now takes `[periodic_u, periodic_v]`. Existing
u-periodic inputs use `[true, false]`; nonperiodic inputs use `[false, false]`.
Both parameter kinds and analytic-section joins retain these row proofs.
A nearly closed row is not a periodic identification.

The [atlas-pair regression](../crates/catallax-shapes/tests/atlas_pair.rs)
cuts a rational unit cylinder by the plane `z = 1/5`. The failing-first
per-chart construction returned eight open branches and no closed branch,
in 545,834 operation units. The joined search returns one closed circle.
It also keeps one open parameter chain when the identical world-space end
rows have not been identified, and rejects an endpoint gap of `1/10^6`.

These operations include chart admission and independent endpoint equation
enclosures. Cells, certificates, endpoint enclosures and work agree at
1, 2, 4 and all available threads. The allowance is unchanged at 20M.

| Periodic cylinder parameter | Pair order | Work |
| --- | --- | ---: |
| u | cylinder, plane | 1,030,454 |
| u | plane, cylinder | 690,222 |
| v | cylinder, plane | 1,755,131 |
| v | plane, cylinder | 1,709,682 |

This is a library proof result. It adds no Boolean-queryable patch-pair curve,
analytic whole-support window certificate or STEP export for these branches.

## Limits

The earlier pair and offset measurements use unit or constant weights.
The joined-cylinder cases above exercise varying positive rational weights.
The poses move both patches together, as R13's did.
A branch along a cell or patch border, and a contact on a border, refuse
rather than answer; certifying them needs boundary-overlap handling, which
does not exist yet. Booleans and STEP writing with spline-with-spline
branches remain open.


Raw evidence retained privately.
