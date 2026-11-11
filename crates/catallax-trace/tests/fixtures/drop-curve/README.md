# Curve drop gate fixtures

`cases.txt` holds outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for `SmSurface::DropAndTrimCurve`
(the drop that SMLib's range answers use for a curve lying in a surface):
the status, each parameter curve as stored (degree, counts, every knot and
the homogeneous control points, `z` and the weight included), the largest
drop distance and the largest normal-line deviation. The library build is
the one described in [`../foundation/README.md`](../foundation/README.md).
Catallax never compiles or links this C++ code.

The cases drop parameter lines of a bicubic bump (along u and v, at a knot,
on its sides, over part of the curve, over a smaller domain the line
leaves, at tolerances of 1e-3, 1e-6 and the default of a thousandth of the
surface's size, with and without keeping far drops); the hill's sides and a
middle line on the hill and on the level plane through its base (the
intersector's range answers), whole and over part of the curve; parameter
lines across and along the C0 knot line of a creased patch; straight and
curved cubics in a level plane, a cubic 1e-4 above it at two tolerances,
and a cubic in a tilted plane; and parameter lines of a collapsed patch
from its pole.

Seven cases are `drefuse`: the port refuses with `SolverRefusal::NotPorted`
and SMLib's answer stays in the raw evidence. Three are side lines outside
a smaller domain, where SMLib finds no start and tries its general curve
drop (`SmSurface::DropCurve`, not ported); one is a cubic 0.1 off the tilted
plane at a control point, which SMLib's plane drop hands to the general
drop; one is a parameter line from the pole and one the pole as a point
curve (a start on a singularity); one is a curve not known to lie on the
surface (the projected side search). The port also refuses a seam start on
a closed surface, the tracer's fallback when the general drop is needed,
and a plane's side lines that meet the curve's box.

One answer deliberately differs from SMLib. `DropIsoCurve` accepts a
parameter line from samples at the curve's knots only, so a curved curve
whose start tangent and knot samples lie on a parameter line (a single-span
quartic in a plane, say) is dropped as that straight line. The port keeps
the line only when a bound over the whole span proves it: the difference
between the curve and the line's image is a spline whose control points
bound it, the line's image being affine on a plane (its values at the
Greville abscissae) and the parameter curve itself on a B-spline surface
(compared when it shares the curve's degree and knots). An unproven line
sends the curve down the general path, which on a plane is its control-point
projection. Every gate case's line is proven; the unit test
`review::a_curved_curve_is_not_dropped_as_a_parameter_line` in
`tests/drop_curve.rs` covers the difference (it failed before the bound).

Binary64 values match within 16 units in the last place. When the fixture
was made every value matched bit for bit. Changing the sample spacing of
`DropIsoCurve`, the keep filter of `DropAndTrimCurve`, the extension of a
parameter line across later breaks or the absolute drop distance changes
cases. These survive: a five-degree angle test for a parameter line, a
tenfold deviation test along it (the cases lie on their lines or far off),
and skipping the trim of a straight drop onto a plane's domain (no straight
case leaves the plane's domain).

To regenerate, build the driver as in the foundation README with
`crates/catallax-trace/tests/fixtures/drop-curve/oracle.cpp`, then run

```sh
python3 crates/catallax-trace/tests/fixtures/drop-curve/generate.py "$ORACLE" "$RAW_DIR"
```

Raw evidence retained privately.
