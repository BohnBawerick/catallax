# Curve-curve gate fixtures

`cases.txt` holds outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for the global curve-curve
intersector `SmCurve::GlobalCurveIntersect` (every solution's type, both
parameters and gap, and for a range its end) with and without its
whole-curve coincidence check, and for that check alone
(`SmCurve::GlobalCoincidenceChecker`, with and without snapping to the
ends: whether more meetings may remain, then its solutions). The library
build is the one described in [`../foundation/README.md`](../foundation/README.md).
Catallax never compiles or links this C++ code.

The cases cross a snake with straight cubics at three points, two of them
its ends, over whole and partial intervals; a wave and a rational wave with
curves through them; two-span and kinked curves (a C0 knot) with a line; a
spiral with its axis; an arc and the closed loop with lines; curves that
miss by 1, 1e-4 and 1e-6. Tangent cases put a parabola on a line through
its vertex, 1e-7 and 1e-4 off it, and an arc on its tangent, at two
tolerances; flat parabolas, a flat cubic and lines at small angles on a line
at tolerances that hold them over a stretch; a seven-span quadratic
B-spline parabola, whole and from its vertex, on a polyline along its
tangent. Touching cases meet end to end, end to interior, 1e-7 off, at C0
knots, and quadratic pieces end to end with parallel tangents. Coincident
cases take a copy of the wave, its reverse and its middle, overlapping
collinear cubics, copies 1e-7 and 1e-3 away, the loop against its copy and
a spiral against its middle, each with and without the coincidence check
and through the check alone; polylines and cubics sharing a straight middle
and parting at both ends; polylines sharing two stretches with a bump of
0.3, 1e-3 and 1e-6 between them; and curves sharing a stretch that covers
neither. Singular cases are a point curve on, off and at the end of a curve
and a cusp crossed at and beside its cusp. A fixed pseudo-random sweep meets
40 cubic pairs, mostly in one plane, at tolerances of 1e-5 and 1e-2.

Five cases are `grefuse`: the port refuses with `SolverRefusal::NotPorted`
and SMLib's answer stays in the raw evidence. Four are pairs of analytic
curves (two pairs of two-point lines, two pairs of arcs) that SMLib meets
with `AnalyticCurveIntersect`; in the fifth the search drops a point onto an
arc, which SMLib hands to the analytic circle (`DropPointFast`), as the
existing curve point drop refuses. The port also refuses a curve
intersected with itself.

Binary64 values match within 16 units in the last place. When the fixture
was made every value matched bit for bit. Dropping the closest-approach keep
of `sm_FilterSolution`, two of the five guesses of a span pair, the
snapping choice of the coincidence check, the extra tangent guess, the merge
of nearby point answers or the early return when the coincidence check
settles the curves changes cases. These survive, unexercised by the cases:
combining separate ranges by their overlap alone, the shift limit of
`sm_RefineCCSolutions`, the first step of
`ClassifyIntersectingRangeSolution`, the discontinuity-only test at span
ends, the closeness of points culled beside a range, and the sample count
of the distance test. No case turns a short or tangent range into a point,
splits a culled zero-gap range into two points, or reaches the end tests of
two discontinuous quadratic spans or the degenerate-curve shortcut.

To regenerate, build the driver as in the foundation README with
`crates/catallax-trace/tests/fixtures/curve-curve/oracle.cpp`, then run

```sh
python3 crates/catallax-trace/tests/fixtures/curve-curve/generate.py "$ORACLE" "$RAW_DIR"
```

Raw evidence retained privately.
