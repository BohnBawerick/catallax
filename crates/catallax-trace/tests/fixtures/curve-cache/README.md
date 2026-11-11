# Curve cache gate fixtures

`cases.txt` holds, for each test curve, the `SmCurveCache` that the pinned
NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, builds through its cache manager
(`SmCurve::CacheMakeOrValidate`). The library build is the one described in
[`../foundation/README.md`](../foundation/README.md). Catallax never compiles or
links this C++ code.

Each line lists SMLib's approximation tolerance, the continuity class at every
unique knot, every leaf span's data and pseudo box, every internal node's
propagated data, and the tree in preorder with each node's box and the span or
internal-node data it refers to. Binary64 values match within 16 units in the
last place because span measures use the arc cosine; every one matched bit for
bit when the fixture was made.

The curves are single Bezier spans of degrees 1, 2, 3 and 5, a straight cubic,
B-splines with smooth, tangent-continuous and kinked joins, non-uniform knots,
a rational full circle with double knots and a weighted cubic. The singular
cases are a curve whose first span collapses to a point and one whose end
tangent vanishes. The near-singular cases are a span 1e-7 long and a curve 3e-9
long. For B-spline curves each Bezier span is the curve itself, so SMLib's
sampled span-to-curve distance is always zero and its tolerance branch is not
reached by this gate.

To regenerate, build the driver as in the foundation README with
`crates/catallax-trace/tests/fixtures/curve-cache/oracle.cpp`, then run

```sh
python3 crates/catallax-trace/tests/fixtures/curve-cache/generate.py "$ORACLE" "$RAW_DIR"
```

Each case runs in a fresh process under a 2 GiB address-space cap, a 20 s CPU
cap and a 40 s timeout. Raw evidence retained privately.
