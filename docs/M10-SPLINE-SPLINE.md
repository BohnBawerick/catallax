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
cell, that exclusion has not removed. Exclusion is exact. It holds when the
boxes round the two control nets are disjoint, or their slabs along either
cell's normal direction are disjoint. It also holds when the two boxes touch
exactly on a side that one patch reaches only on its cell border.

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
point. A cell still unsettled after 48 splits refuses with
`intersect.spline-pair-unresolved` and names the cell and its cover. No
near-enough answer is ever accepted.

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
| M25 wavy pair | 2 open branches | 63 cells, 17 arcs, depth 5 | 3,493,279 | 3,552,719 | 3,503,776 | 3,607,337 |
| Loop, `c = 3/10` | 1 loop | 351 cells, 88 arcs, depth 8 | 16,459,340 | 16,523,717 | 16,466,560 | 16,488,444 |
| Contact, `c = 225/256` | Tangent refusal | refused after 2 splits | 902,655 | 905,946 | 909,875 | 902,655 |
| Small loop, `c = 225/256 - 1/10^6` | 1 loop | 363 cells, 51 arcs, depth 16 | 16,625,889 | 16,711,392 | 16,633,109 | 16,793,100 |
| Disjoint, `c = 225/256 + 1/10^6` | Proved empty | 147 cells, 0 arcs, depth 12 | 3,439,528 | 3,469,471 | 3,446,748 | 3,475,426 |
| Open control | 1 open branch | 12 cells, 3 arcs, depth 1 | 1,020,820 | 1,026,399 | 1,028,040 | 1,021,702 |
| Tiny loop, `c = 225/256 - 1/10^10` | 1 loop | 372 cells, 26 arcs, depth 24 | 16,312,975 | 16,403,746 | 16,320,195 | 16,516,532 |

Cells are the first patch's cells examined, and arcs are certified arcs.
Every row stays inside the allowance; the largest is the scaled small loop at
16,793,100 units. The contact refusal names the exact hilltop at parameters
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

## Limits

No input here has varying weights. The homogeneous formulas cover them,
and the constant-weight copy runs the code path, but no truly rational
surface is measured. The poses move both patches together, as R13's did. Booleans and STEP writing with spline-with-spline
branches remain open.

Raw evidence retained privately.
