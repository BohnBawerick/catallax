# Global point-curve gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for `SmCurve::GlobalPointSolve`
(every solution's type, parameter and value) and `SmCurve::DropPoint` (its
success flag, parameter and distance) on B-spline curves. The library build is
the one described in [`../foundation/README.md`](../foundation/README.md).
Catallax never compiles or links this C++ code.

The curves are a single-span cubic, its rational copy, a three-span spiral, a
closed loop whose ends meet with matching tangents, a closed loop with a
corner at its ends, a two-segment polyline, a parabola, a straight cubic, a
cubic collapsed to a point, a two-point line, a rational circle, and the
arcs listed below. The cases
find nearest, farthest, normal (all and single) and meeting points from points
on the curves, at and beyond their ends, at corners, on the axes of the loops
and just off the smooth loop's seam, over whole and partial intervals; the
circle's axis, where every point is equidistant and SMLib answers with a
range; a point on a line segment; and a fixed pseudo-random sweep of 40
points. The drops run from good
and poor guesses, without a guess (the global solve), normalising, meeting
within tolerances of 1e-5 and 0.5, and snap to either end of the cubic and the
loop with the inward direction along and against the curve.

SMLib minimises and meets on a two-point line with its analytic `SmLine`
(`SmBSplineCurve::CreateAnalyticCurve`, `SmLine::DropPointFast`): the foot of
the point, clamped to the interval, then refined on the spline when the
point is off the line. Those cases put points on, beside, beyond either end
of and on the extension of a segment, over the whole and sub-intervals that
clamp the foot, at meeting tolerances either side of the distance, on a
slanted line far from the origin over the parameter range [2, 5], on a line
shorter than the scaled zero (searched as a spline) and on a rational line
whose unequal weights fail the line test (also searched as a spline). The
point drops without a guess reach the same solve. Removing the clamp, the
refinement or its improvement test changes cases; the end check of the
analytic line is unobserved (a polynomial two-point line always passes it).

SMLib minimises and meets on a circular arc with its analytic `SmCircle`
(`SmBSplineCurve::CreateAnalyticCurve`, `SmEllipse::DropPointFast`). It takes
the point's angle about the circle's axis, clamps it to the arc, and converts it
to the parameter of the spline with doubled interior knots through the polar
converter. Where the conversion is not exact, it refines that parameter by
point inversion. The answer is then refined on the spline when the point is off
the circle. A point on the axis answers both ends of the arc. A point near the
seam of a whole circle answers both sides.

The arc cases use these curves:
- the rational circle;
- a quarter arc of radius 2 in a tilted plane;
- a half arc over the parameter range [2, 5];
- a three-quarter arc with uneven knots;
- a 120-degree arc reparameterised by a Moebius change of weights, so no span converts exactly;
- a 120-degree one-span arc with a knot inserted once, which SMLib doubles;
- a quarter ellipse, which fails the arc test and is searched as a spline.

The test points lie on the arcs, inside and outside them, beyond either end
and on the axis. Near the circle's seam they sit at 1e-7 on either side and
1e-6 outside. The cases include meeting tolerances either side of a 1e-6 gap,
single and all answers, and three sub-intervals of the circle and two
sub-intervals of the other arcs. Drops without a guess reach the same solve.

Sending every arc through the spline search instead changes 27 cases. These
mutations each change cases:
- skipping the knot doubling;
- one angle at the seam;
- one end for a point on the axis;
- always keeping the first of two single answers;
- skipping the point inversion;
- a linear angle guess;
- no clamp to the search interval;
- no meeting filter.

These survive, and the gate does not observe them:
- the period shift of the polar converter's angles;
- the start snap of a closed circle's end in `STEPInversion`;
- the radius test of `IsArc` (without it the quarter ellipse is taken as a circle, and its answers do not change);
- the end check of the circle against the arc (the circle's spline is the arc with its knots doubled, so it always passes).

The arc branch for a cache interval narrower than the curve is not reached,
because `GlobalPointSolve`'s cache spans the whole curve.

The rational two-point line now runs its sampled `SmLine` branch. A circle
whose doubled spline the polar converter cannot take still returns
`SolverRefusal::NotPorted`, because SMLib rebuilds it with `MakeNurb`.
None of these fixture calls reaches that rebuild branch.

The expanded fixture has 1,733 cases. Of its 2,696 binary64 outputs, 2,690
match bit for bit and six differ by one ULP, within the unchanged 16-ULP
candidate allowance. Those six values on the half arc come from `libm`'s
trigonometry. Candidate bits are identical with 1, 2 and 4 workers.
Removing the classification of a closed
curve's seam, the snap direction test, the first-guess exit, the repeated
answer test, the coincident-guess collection, the range assembly or the
boundary angle test changes cases. The bounding-box prune, the line test of
the local solve and the cone cull of normalising only skip work: removing them
leaves every output unchanged.

To regenerate, build the driver as in the foundation README with
`crates/catallax-trace/tests/fixtures/curve-global/oracle.cpp`, then run

```sh
python3 crates/catallax-trace/tests/fixtures/curve-global/generate.py "$ORACLE" "$RAW_DIR"
```

Raw evidence retained privately.

## Rational two-point lines

The analytic factory applies the same sampled line test to rational and
polynomial degree-one curves. Equal weights produce the same bounded
analytic line. Unequal weights that fail the constant-speed test retain the
general spline search. A rational flag alone is no longer a refusal.
The added cases vary six weight scales, four weight ratios, three natural
intervals, partial intervals, endpoint clamps, minimize/intersect requests
and point drops. All results remain unproved candidates.

