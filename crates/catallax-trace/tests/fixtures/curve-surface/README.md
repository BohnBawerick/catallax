# Curve-surface gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for the local curve-surface
solvers `SmSurface::LocalCurveIntersect` and `SmSurface::LocalCurveSolve`
(minimising, maximising, normalising), for `SmSurface::IsSingularity` with
`GetSingularities`, for the curve-surface `SmSurface::SimpleCoincidenceChecker`, and for the
global intersector `SmSurface::GlobalCurveIntersect` (every solution's type,
parameters and deviation), on B-spline curves and surfaces. The library build is the
one described in [`../foundation/README.md`](../foundation/README.md). Catallax
never compiles or links this C++ code.

The surfaces are those of the surface point solver gate. The cases intersect
lines and a cubic with the bump surfaces from good and poor guesses, from
beyond a side and from above, find nearest, farthest and normal pairs, meet
the cylinder with lines tangent to it, 1e-8 outside and inside it, along a
ruling and across it, run lines into the pole of the collapsed patch (where
SMLib moves the answer onto the pole's side of the domain) and test
singularity at, 1e-9 and 1e-4 from, and away from the pole, with and without
the pole-point test. The coincidence cases test a ruling, a full circle and an
arc of the cylinder, an arc 1e-6 off it against tolerances of 1e-5 and 1e-7,
a line half on it (found only when partial coincidence is allowed), a crossing
line and a cubic across a bump, each with and without partial coincidence.

The global cases intersect cubics with the bumps over whole and partial
domains and intervals, a closed loop and a spiral with the cylinder, a cubic
crossing the cylinder twice, a slightly bent curve 0.999 and 0.99999 from its
axis, curves along a knot line of the bump, across its edge exactly at an edge
point and across its corner, point curves on and off a surface, and a fixed
pseudo-random sweep of 48 cubics against the bumps, the loft and the cylinder.
The curve through the pole reaches the singular point inversion, which drops
the point on the iso-curve through the pole. The case listed as `grefuse`
must refuse with `SolverRefusal::NotPorted`: an arc 1e-4 off the cylinder that
SMLib answers with a range, whose classification needs the surface point drop
(`SmSurface::DropPoint`). SMLib's answer for it is in the raw evidence.

Straight curves go to `SmSurface::GlobalLineIntersect`, the line search over
the surface's patches (`SmLSIGlobalSolver`), which these cases now reach:
lines across, through and beside the bumps (over sub-intervals, a parameter
range `[2, 5]`, a sub-domain and as a collinear cubic), across the cylinder,
along its ruling, tangent to it and down its axis, across the loft and the
pole patch and down the pole, and lines placed from the bump's own surface
points (vertical through a point on the knot line u = 1, ending on the
surface, and skimming 1e-3 and 1e-2 under its height near its top). A line
1e-6 inside the cylinder's ruling gets ranges: neighbouring meetings the line
joins along the surface. Curves lying on the bump (its sides and interior
parameter lines, over whole and partial intervals, at two tolerances) pass
SMLib's coincidence check and are answered with ranges built by
`DropAndTrimCurve`. Mutating the polar-box test that stops a patch's
search, the solution merge, the C0 sides, the guess retries or the bound's
expansion changes cases. These survive, unobserved: the bounding-sphere cull
(it only saves work), the duplicate merge's parameter tolerance, the search
from the line's own ends (Newton steps find those points first) and the
second guess from classifying a meeting.

Two answers of the global intersector deliberately differ from SMLib, and the
unit tests in `tests/curve_surface.rs` cover them. When a degenerate curve's
surface point search fails, SMLib ignores the status and keeps what it found
(for a constant curve at the origin on a rational net collapsed there, the
point `(0, 0, 0)`); the port returns the typed failure. The end-point fallback
of the local solve searches the natural domain, and SMLib keeps a point outside
the requested domain (a curve beside the requested half of the unit plane gets
`(0, 0.75, 0.5)`); the port drops it. The point-curve answers also declare three
variables `t, u, v` where SMLib leaves two; the gate does not compare that count.

Binary64 values match within 16 units in the last place; every one matched bit
for bit when the fixture was made.

To regenerate, build the driver as in the foundation README with
`crates/catallax-trace/tests/fixtures/curve-surface/oracle.cpp`, then run

```sh
python3 crates/catallax-trace/tests/fixtures/curve-surface/generate.py "$ORACLE" "$RAW_DIR"
```

Raw evidence retained privately.
