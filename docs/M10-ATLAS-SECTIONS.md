# M10 atlas sections and meetings

The second slice's Boolean target is a box against a lofted body. A smooth
loft face is a `RationalAtlas`: joined exact charts in the source surface's
own parameters. Until now every atlas pair refused with
`intersect.unsupported-pair`. This step adds the three pieces the M25 scope
study listed for that Boolean: atlas sections (its item 5), branch identity
(part of item 3) and meetings with a third surface (item 4). The Boolean now
passes its intersection stage and refuses at two later stages, named below.

## Atlas sections

`RationalAtlas::sections` runs the complete cell search of
[the plane slice](M10-PLANE-SECTIONS.md) on every chart, with the cut
polynomial of [the cylinder and sphere slice](M10-QUADRIC-SECTIONS.md). Arc
ends join by exact equality of their shared parameters, as span joins already
did inside one patch. A periodic face identifies its `u` end rows. An open end
must lie on a chart side with no joined neighbour, which is a physical face
border. Each component is one `PatchSectionBranch` that keeps every chart and
indexes its cells by chart and span. A cut containing a whole chart boundary
row refuses with `intersect.atlas-section-row` before searching.

On the plane slice's bicubic hill, partitioned into one atlas of nine charts,
every plane, cylinder and sphere cut gives the same components, cell counts
and depths as the single-patch search. Tangent contacts refuse with the same
codes and touch points.

## Branch identity

The Boolean asks whether two curves are the same point set.
`PatchSectionBranch::same_component` answers without comparing 3D points:

- Equal chart and surface definitions: the search is exact and deterministic,
  so one component has the same cell ends in both. One shared cell end decides.
- Equal charts, different cutting surfaces: a cell end of one branch where the
  other's cut polynomial is nonzero proves the sets differ.
- A branch and a line or B-spline: a cell-end witness off the other curve
  proves them different.

Anything else still refuses with `intersect.parametric-coincidence`.

## Meetings with a third surface

`PatchSectionBranch::surface_hits` returns the exact transverse meetings with
a plane, cylinder or sphere in traversal order, or `SectionHits::Contained`.

- Meetings at cell ends, including chart seams, are decided by the exact sign
  of the third cut polynomial at those ends.
- Inside each arc cell the two cut polynomials form a square rational system.
  A box where either has a single-signed Bernstein net holds no meeting. Other
  boxes split to depth 10, then try a Krawczyk certificate on the box grown by
  a quarter of its width, so a root on its edge is interior. Certified root
  boxes strictly inside the cell are new meetings. A box touching a known
  meeting is certified together with it; one root there is the known one. A
  known meeting may stand in only if it solves this cell's two equations:
  one from the same chart span does, and one from another span must be
  proved a root first by exact signs. Otherwise a distinct meeting nearby in
  this span would be discarded for it. No equality of isolated roots is ever
  needed.
- A new meeting's coordinates are exact univariate algebraic numbers: the real
  roots of each parameter's eliminant (the resultant of the two equations),
  one inside the certified box on each axis. Points built from isolated system
  roots instead were expensive to enclose for every later consumer.
- Proportional cut polynomials over every cell mean the branch lies on the
  surface: `Contained`. Over some cells only, or a tangent meeting that never
  certifies by depth 24, refuses with `intersect.patch-section-hit`.

Each Krawczyk attempt runs at 128 bits under its own 100,000-unit decision
cap. The numbers crate charges a certificate iteration by term count, not by
rational size. On these systems an uncapped attempt could shrink its box in
tiny steps and spend most of the 20M allowance, so the cap makes such a box
inconclusive and it splits. The request's operation allowance still pays
every unit.

## Verification

```sh
cargo test -p catallax-intersect --test patch_atlas -- --nocapture
cargo test -p catallax-features --test m10_box_loft -- --nocapture
```

Sections of the hill as one atlas of nine charts, including atlas admission.
Every component spans all nine charts. The cell counts and depths equal the
single-patch rows of the earlier slices. The allowance remains 20,000,000
units for every row.

| Cut | Result | Operation units |
| --- | --- | ---: |
| Plane `z = 3/10` | 1 closed loop, 10 arc and 15 empty cells, depth 2 | 213,468 |
| Plane 1e-6 below the hilltop | 1 closed loop, 12 arc and 76 empty cells, depth 21 | 235,875 |
| Plane at the hilltop | `intersect.patch-plane-tangent` naming (5/2, 5/2, 225/256) | 123,198 |
| Cylinder, lowest point 0.3 | 1 closed loop, 10 arc and 5 empty cells, depth 2 | 1,163,542 |
| Cylinder 1e-6 below the hilltop | 1 closed loop, 12 arc and 80 empty cells, depth 22 | 989,425 |
| Cylinder at the hilltop | `intersect.patch-cylinder-tangent`, same point | 271,366 |
| Sphere, lowest point 0.3 | 1 closed loop, 12 arc and 10 empty cells, depth 3 | 1,025,102 |
| Sphere 1e-6 below the hilltop | 1 closed loop, 8 arc and 84 empty cells, depth 22 | 951,006 |
| Sphere at the hilltop | `intersect.patch-sphere-tangent`, same point | 329,011 |

Pair queries on the atlas (orientation in both input orders, membership,
cyclic order, tangents and branch identity) stay under 1,423,000 units per
case. The 1-, 2- and 4-worker runs agree on sections, refusals, meetings and
work units.

Meetings of the hill's `z = 3/10` loop, which crosses chart seams. The
`x = 23/12` plane meets it exactly on the seam `u = 1`. All rows include atlas
admission and the section search, inside the unchanged 20M allowance.

| Third surface | Meetings | Operation units |
| --- | --- | ---: |
| Plane `x = 5/2` | 2, at `y` = 1.2495 and 3.7505 | 775,056 |
| Plane `y = 5/2` | 2, at `x` = 1.2495 and 3.7505 | 932,492 |
| Plane `x = 23/12`, on a chart seam | 2 | 752,620 |
| Plane `z = 1`, above the loop | 0 | 224,814 |
| Sphere touching `z = 3/10` inside the loop | 0 | 2,139,202 |
| The cutting plane with the opposite normal | Contained | |

Two charts joined C1 at `u = 0`, with `z = u - e + 2u²/e` on the left and
`z = u - e` on the right for `e = 1/1024`, meet `z = 0` along `y = 1/10` at
`u = -e` and `u = +e`, one in each span. The left polynomial's other root,
`+e/2`, lies outside its span. Both meetings are kept in traversal order for
the branch and for its reverse, at 665,592 and 665,646 units.

The loft in `m10_box_loft.rs` has a quadratic roof over the floor from (4,0)
to (0,0), peaking at y = 2, 3 and 2 at heights 0, 1 and 2. Its roof is one
periodic atlas face of four charts. Box planes cut it, and the branches meet
the other box planes. Rows include building the loft.

| Cut | Third plane | Meetings per branch | Operation units |
| --- | --- | --- | ---: |
| `y = 1`, two open branches | `z = 1/2` | 1 and 1, at `x` = 0.4152 and 3.5848 | 159,761 |
| `y = 1` | `x = 1` | none | 119,200 |
| `y = 5/2`, one closed loop over the ridge | `z = 1/2` | 2 | 202,813 |
| `y = 5/2` | `x = 2`, on the chart seam | 2, at `z` = 0.3473 and 1.6527 | 261,603 |
| `z = 1/2`, one open branch | `y = 1` | 2, the same two points as the first row | 103,457 |
| `z = 1/2` | `x = 3` | 1 | 71,410 |

The test checks every meeting against both planes through 64-bit enclosures.
The first and fifth rows find the same points from opposite sides.

## The Boolean of a box with the loft

At the atlas-section landing, each case builds the loft and a block through
the recipe evaluator and asks for the Boolean. All six refuse inside the
unchanged 20M allowance. The 1-, 2- and 4-worker runs agree on refusal
details and work units.

| Box | Operations | First refusal | Operation units |
| --- | --- | --- | ---: |
| 2 x 10 x 1 at (1, 1, 1/2): crosses the roof's sides | subtract, union, intersect | `boolean.imprint-ambiguous`: curve does not provide a supported elimination chart | 1,074,776 |
| 10 x 10 x 1 at (-3, 5/2, 1/2): takes the top off the ridge | subtract, union, intersect | `shapes.parametric-query`: a joined atlas needs a trimmed body side query | 2,447,617 |

Before this step all six refused at once with `intersect.unsupported-pair`.
The intersection stage now finishes: atlas sections, branch identity, the
meetings and exact identity of their shared vertices all hold. A work profile
of an earlier version, whose meeting points were isolated system roots, spent
18.1M and 19.1M units in one exact zero proof that two such vertices
coincide. With eliminant coordinates the side-box case's eight such proofs
cost 26,836 units together.

The next missing pieces are therefore named by the refusals:

1. An imprint chart for a patch branch. Imprinting a curve on a face asks the
   curve for an elimination chart; a branch can supply one through its patch
   parameters and cells.
2. A side query against a trimmed atlas face, to classify points against the
   lofted body.

The later [imprint query slice](M10-IMPRINT-QUERIES.md) moves these refusals
without yet proving a Boolean solid.

## Limits

Tangent meetings refuse. There is no STEP writer for patch branches. The
candidate tube certificate still covers open single-span branches only, so a
branch crossing chart seams keeps its exact definition without a checked
representation. [Cone and ring-torus equations](M10-CONE-TORUS-SECTIONS.md)
extend the shared search and meeting checks. Spline against spline remains open.
