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

Further operations exercise the tracing stages at chosen parameter pairs:
`ComputePointValues` alone and after a previous point (`point`), one
`ComputeStepSize` and `ComputeNextPoint` step from the meeting point near a
pair (`step`), a whole `TraceCurve` with no start, end or through points
(`trace`), and `ClassifyBoundaryStartPointIntersection` at meeting points on
the sides (`classify`).

The `point` cases include parallel normals where one surface is planar
(`ComputeSingularityPoint`): a wall along the cylinder patch's middle ruling,
facing either way, against the patch and its copy with v reversed (normal
inward), alone, with either surface first, after a point on the ruling and
after points off the meeting whose tangents point either way (a tangent
curve), and the hill's plateau on its level (a tangent point).

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

Six `boundary` cases are `brefuse`: the port refuses with
`SolverRefusal::NotPorted`, and the driver records `-1` with SMLib's answer
kept in the raw evidence. They are the loft against its sloping plane at both
tolerances (its u = 0 side is straight, which SMLib intersects with
`GlobalLineIntersect`); the hill on the plane through its base (its sides lie in the
plane, giving range answers); the bump against a plane domain that ends inside
the meeting box (the plane's side lines, `SmLine`, would be intersected); and
the patch with a pole against both planes (its last side is straight, which
SMLib intersects with `GlobalLineIntersect`).

The `intersect` cases add walls 1e-3, 1e-2 and 1e-1 inside the cylinder
patch's ruling, levels 1e-2, 3e-2 and 1e-1 below the bump's highest point,
and curved sheets that meet only at a corner, from outside and from inside
the corner's sector. Eleven of them are `irefuse`, the driver recording `-1`
with SMLib's curves kept in the raw evidence: the loft, the patch with a pole
and the flat patches have a straight side (`GlobalLineIntersect`); the hill's
sides lie in the plane through its base (range answers); the bounded plane's
side lines cross the meeting box (`SmLine`). Three cases run no `intersect`:
SMLib's `DoIntersection` crashes (a segmentation fault) on the flat patches
that share an edge, lie 1e-9 apart, or share half an edge.

Binary64 values match within 16 units in the last place. When the fixture was
made every value matched bit for bit except five tangent-plane angles (three
on traced points, two on side start points), one unit apart: the angle is an
`acos`, and the port's `libm` and the C library round differently. Swapping the subdivision order of
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
a left-handed pair). In the whole intersection, dropping the merge of a short
last span, the two-degree test on the normals at interior meetings, the
four-tolerance test on touching curves, the
mid-span test of `IsPointOnCurve`, the order of the interior guesses, the
interior search itself or the grazing-box process changes cases. These
survive, unexercised by the cases: the near-tangent curve type (no crossing
curve with more than one span stays within five degrees), the trace behind a
start that can run both ways, the five-degree interior test, the
one-tolerance test of `IsPointOnCurve`, the removal of a touching curve lying
on another, the corner-sector test, the second-point check of a trace against
earlier curves, the move of the larger patch's midpoint onto its tangent
plane, the uniqueness test for a one-point curve, and the through point kept
for a start that gave no curve.

To regenerate, build the driver as in the foundation README with
`crates/catallax-trace/tests/fixtures/intersector/oracle.cpp`, then run

```sh
python3 crates/catallax-trace/tests/fixtures/intersector/generate.py "$ORACLE" "$RAW_DIR"
```

Raw evidence retained privately.
