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

Two arc-pair cases are `grefuse`: the port refuses with
`SolverRefusal::NotPorted` and SMLib's answer stays in the raw evidence.
The two formerly refused polynomial two-point line pairs now pass through
`AnalyticCurveIntersect`'s line branch (`SmLine::IntersectWithLine`). A search that drops a point onto an arc goes to
the analytic circle (`DropPointFast`) through the curve point drop, as in
SMLib; one case, an arc against a straight cubic, does so. The port also
refuses a curve intersected with itself.

## Analytic line pairs

The expanded fixture has 468 cases, including 288 added line-pair calls.
Sixteen known SMLib misses are excluded from the 452-case correctness count.
Both curve orders and both coincidence-check flags run on each added pair.
The lines cross in the interior and at ends, touch, overlap in both parameter
orientations, coincide, miss in 3D and have small or collapsed requested
intervals. Whole and partial intervals, natural parameter ranges away from
`[0,1]`, gaps on both sides of the tolerance and seeded random 3D crossings
exercise the cache's analytic dispatch.

The segment helper checks endpoint proximity before its interior closest-line
solve, keeps at most two distinct points and retains SMLib's near-parallel
range rule. Range endpoints sort along the first curve, keeping the second
curve's direction. The cache branch returns those ranges without the extra
range refinement used for explicitly analytic input objects in SMLib.

The determinism check runs the complete fixture on 1, 2 and 4 workers and compares
candidate bits in their original order. The public-API regression tests
check the formerly refused crossing and a reversed coincident range.
All returned points and ranges remain unproved candidates. Exact acceptance
and deterministic work allowances stay unchanged.

Binary64 values match within the unchanged 16-ULP comparison allowance.
All 1,692 values in the expanded fixture match bit for bit. Dropping the closest-approach keep
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

## R1 known misses on line pairs

Eight `known-miss` rows cross a segment from `(0,0,0)` to `(L,0,0)`
with one from `(L/2,-L,0)` to `(L/2,L,0)`, for `L = 1e-10` and `1e-8`,
both curve orders and both coincidence-check flags. The exact dyadic control
points define the same point at `t1 = t2 = 1/2`. Pinned C++ SMLib nevertheless
returns success with zero candidates. Its closest-line helper rejects the
cross product at its absolute `SM_EFF_ZERO_SQ` threshold. The `1e-6` copies
return the expected midpoint candidate.

Eight further rows collapse one requested interval to `[1/2,1/2]` at the
common midpoint of two unit-scale crossing lines. That closed interval
contains the meeting, but SMLib's analytic line branch returns no candidate
when either evaluated segment has zero length. The point witness is checked
in both curve orders and with both coincidence-check flags.

The port faithfully reproduces those misses. They are explicit R1 exceptions,
excluded from the ordinary correctness count. Separate tests check their
classification, candidate parity and common point witnesses. Empty candidate
lists cannot certify an empty intersection; exact proof remains separate.

To regenerate, build the driver as in the foundation README with
`crates/catallax-trace/tests/fixtures/curve-curve/oracle.cpp`, then run

```sh
python3 crates/catallax-trace/tests/fixtures/curve-curve/generate.py "$ORACLE" "$RAW_DIR"
```

Raw evidence retained privately.
