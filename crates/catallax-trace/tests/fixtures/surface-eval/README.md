# Surface evaluation gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for `SmBSplineSurface::Evaluate`,
`EvaluatePoint` and `EvaluateSimple`. The library build is the one described
in [`../foundation/README.md`](../foundation/README.md). Catallax never compiles
or links this C++ code.

The surfaces are a bicubic bump with interior knots (polynomial and rational),
a rational loft-like patch, a cubic patch whose u = 0 edge collapses to a pole,
and a bilinear patch. The cases evaluate points and derivatives up to fourth
order in full and upper-triangular form, at interior knots from each side, at
corners and outside the domain (clamped), and at the pole with SMLib's repair of
vanishing tangents, which re-evaluates at second order and, when that fails,
samples nearby points, and surface normals with their stepping, second-derivative and sampling fallbacks at the pole. Requests beyond fourth order and upper-triangular
requests with unequal orders fail as in SMLib. Point evaluation has no tangent
or touching meeting, so this gate's special cases are the singular pole, its
near-singular neighbourhood, domain boundaries and knots.

To regenerate, build the driver as in the foundation README with
`crates/catallax-trace/tests/fixtures/surface-eval/oracle.cpp`, then run

```sh
python3 crates/catallax-trace/tests/fixtures/surface-eval/generate.py "$ORACLE" "$RAW_DIR"
```

Raw evidence retained privately.
