# Plane gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for `SmPlane`: the bilinear net
`MakeNurb` builds over the analytic domain, `EvaluatePoint`, `Evaluate` and
`EvaluateSimple`, `EvaluateNormal`, `IsSingularity`, `IsClosed`, `IsOnSeam` and
`DropPointFast`, and for `SmSurface::GlobalPointSolve` and
`SmSurface::GlobalCurveIntersect` on planes. The library build is the one
described in [`../foundation/README.md`](../foundation/README.md). Catallax
never compiles or links this C++ code.

The planes are two made by `SmPlane::CreateCanonical` (unbounded, over
`±SM_INFINITE_PARAMETER`; one with unit axes and one with non-unit,
non-perpendicular axes) and four bounded ones with scaled parameters: axes 1e-3 off
perpendicular (made perpendicular by `SetCanonical`), axes 1e-13 off (kept),
power-of-two scales, and scales 0.3 and 1.7 on rotated axes. Evaluation runs at
interior points, corners, points outside the domain (clamped) and 1e-9 from
the origin, for derivative orders up to 3 by 2, full and upper-triangular. The
point drops minimise, normalise and meet from points on, off, 1e-9 off and
beyond the planes, inside and outside sub-domains and a domain larger than the
plane's, with and without a target distance. The global point solves run the
same points. The curve cases cut the unbounded plane with a crossing cubic, a
quadratic tangent to it, one 1e-7 above it, a cubic touching it at an end, a
point curve on it and a two-span cubic crossing it four times, over the whole
domain and two sub-domains, and cut the tilted and a bounded plane.

Straight curves reach `SmPlane::GlobalLineIntersect`: lines crossing the
plane inside and outside the curve's interval and the domain, crossing it at
the curve's end, a collinear cubic, a line over the parameters `[2, 5]`, lines
across the tilted, a bounded and the scaled plane, and lines parallel to the
plane 1 and 1e-4 above it and beside a sub-domain. Lines lying in the plane,
1e-6 and 1e-9 above it, or crossing a sub-domain in it, get a range answer:
both ends of the line's piece over the domain, with their parameters and
distances (each solution's end is printed for ranges). Removing the line's
meeting tolerance, its interval test, the domain clamp, or the parallel
line's domain-overlap test fails cases; removing the parallel line's start
distance test does not, because the end distance test rejects the same lines.

Binary64 values match within 16 units in the last place; every one matched bit
for bit when the fixture was made. Evaluating the point with the scale applied
after the axis (as `Evaluate` does) instead of before it (as `EvaluatePoint`
does) changes 2 values by up to 2 ulp, and skipping `DropPointFast` in the
global point solve changes 54 values by up to 16 ulp; the gate reports both as
lost bit-identity. Removing the domain clamp, the meeting tolerance or the axis
orthogonalisation fails cases. The `IsInDomain` filter after `DropPointFast`
only rejects drops outside the natural domain, and for the plane larger than
its domain the B-spline search that follows finds the same point, so removing
the filter leaves every output unchanged.

The port refuses inputs SMLib accepts without checking, with
`Fault::InvalidInput`: planes with a zero or non-finite parameter scale (whose
projection would divide by zero) or non-finite frame or domain, drops whose
projection or distance is not finite, and derivative orders beyond the
evaluator's limits (as for other surfaces). The unit tests in
`tests/plane.rs` cover these. A drop over a domain larger than the plane's
reads the plane's analytic extension, as SMLib's does, while point
evaluation clamps to the plane's own domain.

One answer of `GlobalLineIntersect` deliberately differs from SMLib. Its
parallel test (0.0001 degrees) admits lines with a slight slope, and SMLib
then returns no meeting as soon as one end of the trimmed line is beyond the
tolerance, even when the line crosses the plane between its ends. The port
answers that crossing when the trimmed ends lie on opposite sides of the
plane, meets nothing only when both ends are beyond the tolerance on one
side (the gap is linear along the line), and refuses with
`SolverRefusal::NotPorted` when only one end is within it. No gate case
reaches the difference; the unit tests in `tests/plane.rs` module
`near_parallel` cover it (the crossing test failed before the change).

To regenerate, build the driver as in the foundation README with
`crates/catallax-trace/tests/fixtures/plane/oracle.cpp`, then run

```sh
python3 crates/catallax-trace/tests/fixtures/plane/generate.py "$ORACLE" "$RAW_DIR"
```

Raw evidence retained privately.
