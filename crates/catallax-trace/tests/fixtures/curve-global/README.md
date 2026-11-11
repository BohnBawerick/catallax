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
cubic collapsed to a point, a two-point line and a rational circle. The cases
find nearest, farthest, normal (all and single) and meeting points from points
on the curves, at and beyond their ends, at corners, on the axes of the loops
and just off the smooth loop's seam, over whole and partial intervals; the
circle's axis, where every point is equidistant and SMLib answers with a
range; a point on a line segment; and a fixed pseudo-random sweep of 40
points. The drops run from good
and poor guesses, without a guess (the global solve), normalising, meeting
within tolerances of 1e-5 and 0.5, and snap to either end of the cubic and the
loop with the inward direction along and against the curve.

Cases listed as `grefuse` and `drefuse` must refuse with
`SolverRefusal::NotPorted`: SMLib solves minimising and meeting on a two-point
line and on a rational quadratic arc with its analytic `SmLine` and `SmCircle`
drops. The port's arc test is conservative: every rational quadratic that is
not degenerate refuses these operations. SMLib's answers for them are in the
raw evidence.

Binary64 values match within 16 units in the last place; every one matched bit
for bit when the fixture was made. Removing the classification of a closed
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
