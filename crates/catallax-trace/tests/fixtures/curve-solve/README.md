# Curve solver gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for `SmCurve::LocalPointSolve`,
`LocalCurveIntersect`, `LocalCurveSolve`, `CurveMaxDistanceBetween`,
`SimpleCoincidenceChecker`, `IsClosed` and `IsPeriodic` on B-spline curves.
The library build is the one described in
[`../foundation/README.md`](../foundation/README.md). Catallax never compiles or
links this C++ code.

The cases drop points onto curves (including the centre of a circle, points far
inside a tight bend, points beyond a curve end, a curve with a collapsed end and
a curve a thousand units long), intersect crossing, tangent and near-tangent
curves (a parabola with lines 1e-8 below and above its vertex), disjoint curves,
curves meeting at an end and overlapping lines, find closest and farthest
pairs, measure the largest gap between a curve and a copy moved by 1e-6, and
check coincidence of overlapping lines, a curve against its own copy on
overlapping intervals, and near and far copies, with and without pseudo boxes.

Binary64 values match within 16 units in the last place; every one matched bit
for bit when the fixture was made. Small changes to SMLib's Jacobian limits
converge to the same roots and do not show in these outputs; removing a limit
changes six cases.

To regenerate, build the driver as in the foundation README with
`crates/catallax-trace/tests/fixtures/curve-solve/oracle.cpp`, then run

```sh
python3 crates/catallax-trace/tests/fixtures/curve-solve/generate.py "$ORACLE" "$RAW_DIR"
```

Raw evidence retained privately.
