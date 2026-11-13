# Curve drop gate fixtures

`cases.txt` holds outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for `SmSurface::DropAndTrimCurve`
(the drop that SMLib's range answers use for a curve lying in a surface,
operations `drop` and `drefuse`) and for `SmSurface::DropCurve` called
directly (`SmPlane::DropCurve` on a plane; operation `general`): the status, each parameter curve as stored (degree, counts, every knot and
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

The 66 `general` cases walk curves onto surfaces with the general drop:
a cubic on the bump and a five-point cubic half a unit over it at
tolerances of 1e-3, 1e-6 and 0.05, with and without keeping far drops, over
parts of the curve and over a smaller domain they leave; a line well above
the bump; curves starting on and running into the bump's sides; a vertical
line over the bump (its tangent along the normal), one five degrees off
vertical and one ending on the bump; a curve 0.1 off the tilted plane and a
curve rising off the level plane (which `SmPlane::DropCurve` hands to the
general drop) and one leaving a small domain of the plane; helices on and
around a rational cylinder across and from its seam, in both parameter
directions, a ruling on its seam and a line normal to it; lines from,
into and past the pole of the collapsed patch with the pole on the `u` and
on the `v` side, and the pole as a point curve; a line across the crease;
lines on a flat bilinear patch, inside, along a side and leaving it;
straight cubics at uneven speed along the flat patch's sides and the
cylinder's seam; a cubic leaving the bump's smaller domain just before its
end; and a domain larger than the surface's (SMLib's error). Over these
cases the instrumented SMLib runs 746 of the 885 code lines of
`DropCurve` and the `SmDropPt` helpers it calls. Not run: the seam-jump
correction (`CheckSeamJump` finds no jump), the dog-leg straightening at a
pole, the rescue of a walk that leaves the domain within the tolerance of
the curve's end, the rescue of a kept drop beyond the tolerance, and the
filter of several curves when far drops are not kept.

Mutations that survive: widening the three-degree tangent test of
`IsGoodTangentDrop` to five degrees, skipping `FixBoundaryTangents` (the
walks along sides already have no cross component), and the clamped
instead of the extrapolated normal in the closing measure (it moves one
deviation of 1.3e-15 by two units in the last place, inside the gate's
tolerance; the port matches SMLib bit for bit with the extrapolated
normal). These fail the gate: two instead of three span samples, no
second-order Taylor guess, no offset term in the dropped tangent, ten
instead of twenty closing samples, and no duplicate removal after the seam
copies.

Two cases remain `drefuse`: the port returns `SolverRefusal::NotPorted`
and SMLib's answer stays in the raw evidence. One is a cubic 0.1 off the
tilted plane at a control point, which reaches the tracer's fallback. The
other is a curve not known to lie on the surface, which needs the projected
side search. Seam starts and the tracer fallback remain typed refusals.

The two formerly refused pole cases now compare complete positive outputs:
a parameter curve leaving the pole and the pole as a point curve. Before
changing the expectations, direct comparison with the retained C++ calls
matched all 28 binary64 outputs and every integer exactly. Singular starts
recover the varying UV coordinate by dropping the nearby curve point. The
separate `../drop-start` gate exercises `ComputePointValues`, including the
nearby tangent-point drop. Plane boundary lines also enter the curve-pair
search. The full fixture is checked with 1, 2 and 4 workers.

One answer deliberately differs from SMLib. `DropIsoCurve` accepts a
parameter line from samples at the curve's knots only, so a curved curve
whose start tangent and knot samples lie on a parameter line (a single-span
quartic in a plane, say) is dropped as that straight line. The port keeps
the line only when a bound over the whole span proves it: the difference
between the curve and the line's image is a spline whose control points
bound it, the line's image being affine on a plane (its values at the
Greville abscissae) and the parameter curve itself on a B-spline surface
(compared when it shares the curve's degree: both are first refined by
knot insertion into their common knot vector, since paired control points
bound a difference only on one basis; ends within the scaled zero are
mapped onto the curve's). An unproven line
sends the curve down the general path, which on a plane is its control-point
projection. Every gate case's line is proven; the unit test
`review::a_curved_curve_is_not_dropped_as_a_parameter_line` in
`tests/drop_curve.rs` covers the difference (it failed before the bound),
and the unit test `close_but_different_knots_do_not_prove_a_parameter_line`
in `surface_drop_curve.rs` covers knots that differ by 2^-34 (it failed
while close knots counted as one basis).

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

After singular-side recovery, all 132 cases pass and all 8,000 binary64
outputs match bit for bit. Only the two documented singular refusal rows
changed to complete positive outputs. Candidate bits match on 1, 2 and 4
workers. The two remaining `drefuse` rows retain typed refusals.
