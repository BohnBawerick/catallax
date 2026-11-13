# Surface intersector gate fixtures

`cases.txt` holds outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for
`SmAdvSurfaceIntersector::DoIntersection` (the intersector that
`SmSurface::GlobalSurfaceIntersect` runs) as a whole and stage by stage, read
through a probe subclass of the intersector: the constructor's tolerance, closed flags and
clipped domains; the root-box test; the `SM_SO_INTERSECTION_TEST` search with
the advanced intersector's overrides (the union of the meeting boxes and
domains of the leaf pairs, and their count); and `DoBoundaryPlaneTest`. The
library build is the one described in
[`../foundation/README.md`](../foundation/README.md). Catallax never compiles or
links this C++ code.

The cases are the eleven surface cases of the M25 scope study's coverage
driver (a bicubic bump and its rational copy against level and tilted planes,
with the default and a tight tolerance; a hill against planes through its
plateau, at its top and 1e-6 below it; a rational cylinder patch against a
wall along its ruling and across it; two bumps; a loft against a sloping plane
with two tolerances), and variations: planes missing the surfaces and touching
them at a corner, an edge and the top, sub-domains, a requested domain beyond
the bump's natural one, a smaller plane domain, swapped surfaces, walls 1e-6
off the cylinder's ruling, the patch with a pole against planes through and
below its pole, more levels, and flat patches side by side with a gap of 5e-4,
2e-2 and 1e-9, sharing an edge, overlapping, at an angle and turned 0.005 and
0.02 radians about z. Each case runs the search alone, with the
boundary-plane test, on through `FindBoundaryStartPoints` (`boundary`: the
start points found on the sides of both surfaces, the through-point count and
the touching curves), and the whole `DoIntersection` (`intersect`: each
curve's type and deviation, then the stored 3D curve and the two parameter
curves, knots and control points as kept). These are Gate B, spline surface
with plane, and spline with spline where the cases have two patches.

The `self` operation runs `DoIntersection` on one surface given twice, which
is how SMLib's `GlobalSurfaceIntersect` meets a surface with itself (the
constructor sees the one pointer, skips the side search and rejects leaf pairs
whose parameter domains overlap). Cases come with the default and a
1e-4 tolerance where the surface meets itself: a planar cubic loop with
control points (0,0), (4,3), (-2,3), (2,0) that crosses itself once, swept
along z (straight, slanted, rational, tilted and bent into a quadratic sweep),
the same loop tilted until its two branches no longer meet (empty results),
domain pairs that share the side at the crossing, are apart or split the loop
at other places, and the bump, its rational copy, the cylinder patch, the loft
and the hill against themselves over the whole domain and over split domains
(no self-intersection, empty results). A quartic loop whose two branches touch
at one point with a common tangent (control points (0,0), (1.847201,3.208322),
(2.464501,y), (-0.007119,2.438416), (1.683092,0), y = 1.2892583905 found to a
2e-15 residual), swept along z, touches itself along a line without crossing:
straight, slanted and on split domains, and with the third control point's y
moved by 1e-6, 1e-3, 0.2 and 0.5 (branches pulled apart or made to cross). In
all, sixty-one cases: twenty-six give one curve and thirty-five give none. SMLib
returns no curve for the exact touch and for the 1e-6 and 1e-3 moves.

Further operations exercise the tracing stages at chosen parameter pairs:
`ComputePointValues` alone and after a previous point (`point`), one
`ComputeStepSize` and `ComputeNextPoint` step from the meeting point near a
pair (`step`), a whole `TraceCurve` with no start, end or through points
(`trace`), and `ClassifyBoundaryStartPointIntersection` at meeting points on
the sides (`classify`).

The `point` and whole-intersection cases also include curved surfaces whose
normals are parallel where they meet: quadratic patches `z = a x^2 + b y^2 +
c x y`, exact as degree-2 nets, meeting at the origin. They are a saddle on a
level plane (four branches), two caps touching (a tangent point), an umbilic
patch against one that bends more along x and the reverse (tangent curves),
opposed saddles, a saddle against a turned patch, patches with v reversed
(opposed normals) and a patch against itself (coincident sides), each at the meeting alone and after points beside it. These
run the Dupin indicatrices of `ComputeSingularityPoint`.

The `point` cases include parallel normals where one surface is planar
(`ComputeSingularityPoint`): a wall along the cylinder patch's middle ruling,
facing either way, against the patch and its copy with v reversed (normal
inward), alone, with either surface first, after a point on the ruling and
after points off the meeting whose tangents point either way (a tangent
curve), and the hill's plateau on its level (a tangent point).

The `step` and `trace` cases include the collapsed patch (its u = 0 side is a
pole) against the wall y = 0 through the pole, with either normal, from
points beside the pole and at the patch's middle: steps and traces that reach
the pole, where `ComputePointValues` with a previous point and a step size
steps back from the singularity (`StepBackFromSingularity`, sixteen times
across the cases). Evaluating the Hermite span at 0.01 instead of 0.99,
taking the other parameter from the stepped-back point, skipping the step
back, keeping the point's own tangent for the span, or keeping the singular
`u` derivative changes cases.

When two side planes graze, SMLib intersects the side curves with the curve
intersector, which is not ported; the port refuses with
`SolverRefusal::NotPorted` there and the driver records `-1`, as for disjoint
root boxes. The `boundary` operation records a grazing leaf pair
(`GrazingBoxProcess`) as `-2` and a grazing side-plane pair as `-3`.

The `boundary` cases add levels at the highest point of the bump's u = 0 side
(the meeting curve runs inside along the side), at the highest and lowest
points of its v = 0 side (the plane touches the domain at one point, a
touching curve) and 1e-9 above and below, levels through the bump's corners,
and two curved patches whose sides lie in the planes y = 0 and y = 3 and cross
there (crossings on both surfaces' sides, found before and after inner ones,
and on a side 1e-11 inside the other's).

Two `boundary` cases meet sides in ranges: the hill on the plane through its
base (its sides lie in the plane) and a quadratic patch against itself
(coincident sides). A long range becomes an intersection curve: the side
trimmed to the range and dropped onto both surfaces with `DropAndTrimCurve`,
typed by `sm_ComputeCurveType` (tangent for the hill, coincident for the
patch); start points on it are removed. Straight sides met with a plane (the
loft's u = 0 side, the pole patch's last side) take the plane's analytic
`GlobalLineIntersect`; a plane's own sides (`SmLine`) and other straight sides
take the line search over the other surface's patches.

The `intersect` cases add walls 1e-3, 1e-2 and 1e-1 inside the cylinder
patch's ruling, levels 1e-2, 3e-2 and 1e-1 below the bump's highest point,
and curved sheets that meet only at a corner, from outside and from inside
the corner's sector. Five of them are `irefuse`, the driver recording `-1`
with SMLib's curves kept in the raw evidence: the flat patches overlap, and
the port refuses before their sides' ranges, at the search stage. Four
leaf pairs lie in a box thinner than twice the tolerance (SMLib's
`GrazingBoxProcess`) and one pair's boundary-plane test is ambiguous; both
are not ported. The general curve drop (`SmSurface::DropCurve`) these
ranges then use is ported and gated in `../drop-curve`. The
hill on its base plane and the patch against itself now match SMLib's whole
intersection, range curves included. Making every curve type tangent where
the normals cross, the umbilic test's coincidence, the gap-ratio filter,
the drop-distance filter, the removal of start points on a range curve or
the long-range threshold survive: the range curves here are tangent or
coincident (from the principal directions), with equal zero gaps at both
ends, and no start point lies on them. Three short ranges take the
short-range branch (the tighter end moved by a local normal solve).
Three cases run no `intersect`:
SMLib's `DoIntersection` crashes (a segmentation fault) on the flat patches
that share an edge, lie 1e-9 apart, or share half an edge.

One classification deliberately differs from SMLib. Its planarity test for
parallel normals (`smsurf_ComputeSurfacePointType`) reads only `N . S_uu`
and `N . S_vv`, so a twisted point such as the saddle `z = u v` passes for
planar, and its meeting with a plane is called a tangent point. The port
also requires `N . S_uv` to vanish; such a point is a saddle, whose meeting
with a plane is a singularity. No gate case reaches the difference; the unit test
`review::a_mixed_curvature_saddle_is_not_planar` in `tests/intersector.rs`
covers both surface orders.

A second differs at a singularity met at a start point with a previous
direction: SMLib measures the branches leaving into both domains against
that direction but, when the first of them is the nearest, keeps the first
branch found, which may leave a domain. The port takes that inside branch.
No gate case reaches the difference; `review::a_start_singularity_takes_an_inside_branch`
covers a saddle met on each side of its domains.

Binary64 values match within 16 units in the last place. When the fixture was
made every value matched bit for bit except eleven tangent-plane angles (five
on traced points, two on stepped points, four on side start points), one unit
apart but for one two units apart: the angle is an `acos`, and the port's
`libm` and the C library round differently; for normals more than a right
angle apart the angle is `pi` less the `acos`, near 2, and one unit there is
two units of the result near 1.14. Swapping the subdivision order of
`sm_FindLargestNonLeaf`, dropping the tolerance expansion of the pseudo-box
test, changing the default tolerance, removing the disjoint outcome of the
plane test or the clipping of a requested domain to the natural one changes
cases. Widening the opposed-normal test of the side planes from 1e-4 to 1e-2
leaves every output unchanged: the turned patches' planes do not separate the
surfaces either way. In `FindBoundaryStartPoints`, adding side crossings
on the other surface's side at once instead of last, dropping the 100 *
EFF_ZERO side tolerance, taking the touching curve's point from the side curve
alone or from the wrong surface, dropping a start point's deviation, using the
domain minimum for a side, or keeping touches within four tolerances of an
earlier one changes cases. Two mutations survive: dropping the tolerance
expansion of a side curve's box (it only admits sides whose answers lie
beyond the tolerance) and keeping touches that lie on a start point (no case
has a touch near a start point). In the planar branch of `ComputeSingularityPoint`,
turning off the alignment with the previous tangent or the reversal for
opposed normals, choosing the other principal direction, dropping the swap of
principal directions by normal curvature, turning tangent points into
refusals, or widening the planar test changes cases; dropping the final
handedness fix of the principal directions survives (no case reaches it with
a left-handed pair). On the curved side, deduplicating singular directions,
the sign of the turn between the indicatrices, the negated candidates after a
previous point and the flatter direction of a coincidence change cases;
dropping the thousandfold clamp of a negligible principal curvature, either
umbilic test, the swap of a saddle's asymptotes for opposed normals, or
widening the inside test of a start point's branches survives. In the whole intersection, dropping the merge of a short
last span, the two-degree test on the normals at interior meetings, the
four-tolerance test on touching curves, the
mid-span test of `IsPointOnCurve`, the order of the interior guesses, the
interior search itself or the grazing-box process changes cases. These
survive, unexercised by the cases: the rejection of leaf pairs on overlapping domains in a self-intersection (removing it changes no case), the near-tangent curve type (no crossing
curve with more than one span stays within five degrees), the trace behind a
start that can run both ways, the five-degree interior test, the
one-tolerance test of `IsPointOnCurve`, the removal of a touching curve lying
on another, the corner-sector test, the second-point check of a trace against
earlier curves, the move of the larger patch's midpoint onto its tangent
plane, the uniqueness test for a one-point curve, and the through point kept
for a start that gave no curve.

## Analytic line-pair follow-up

The polynomial line-pair port removes five stale whole-intersection refusal
expectations at generator case indices 34, 36, 37, 39 and 40. Before changing
the gate, each original request ran directly through the pinned C++ probe
and the port. All 480 binary64 values and all integer outputs matched bit for
bit, including curve types, deviations, 3D curves and both parameter curves.
The candidate counts were 1, 1, 4, 3 and 3 respectively. The gate now compares
the complete outputs instead of requiring `SolverRefusal::NotPorted`.
The three recorded C++ crashes at indices 33, 35 and 38 stay excluded.

These are unproved tracer outputs. In particular, the parallel-patch case
separated by 1e-3 returns tolerance-based candidates; matching those outputs
is not proof that the exact surfaces intersect. Catallax's proof rules and
accuracy caps remain unchanged. The full gate also compares ordered candidate
bits with 1, 2 and 4 workers.

To regenerate, build the driver as in the foundation README with
`crates/catallax-trace/tests/fixtures/intersector/oracle.cpp`, then run

```sh
python3 crates/catallax-trace/tests/fixtures/intersector/generate.py "$ORACLE" "$RAW_DIR"
```

Raw evidence retained privately.
