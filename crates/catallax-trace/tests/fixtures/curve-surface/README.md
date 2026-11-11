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
the point on the iso-curve through the pole. The arc 1e-4 off the cylinder
is listed as `grefuse` and must refuse with `SolverRefusal::NotPorted`: SMLib
answers it with a range through its coincidence check, whose curve drop meets
two analytic curves (`AnalyticCurveIntersect`). SMLib's answer for it is in
the raw evidence.

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

Range answers reach the solution combiner
(`sm_GlobalCurveIntersectCombineSolutions`) through the coincidence checks
along patch sides that are at most C0. Forty cases put a curve on such
a side over a short stretch between the whole-curve checker's samples. The
sides are the v = 0 edge and the crease of two folded bilinear sheets, and
the bump's v = 0 edge, where each stretch is cut exactly from the edge curve
by its blossoms. Off the stretch the curves rise off the surface or run
outside the domain. Together they cover:

- ranges kept as found;
- ranges past the 2% coincidence threshold, rebuilt from `DropAndTrimCurve`
  (with the points where the curve is tangent to the surface dropped), and
  the rebuild abandoned when the curve leaves the surface;
- consecutive ranges whose gap drops to the surface within the tolerance and
  join, and gaps lifted 2e-5, 3e-5 or 1e-3 off it that stay apart;
- ranges 5e-6 above the edge or tilted to 8e-6 off it, which the touching
  classification (`ClassifyIntersectingRangeSolution`) turns into points or
  keeps;
- a range 1e-9 long, and a curve that stalls on the edge.

Turning off the gap join, the touching classification, the rebuild, its
tangent-point filter or the range merge each fails between 1 and 11 cases.
The combiner's point-beside-range cull is unobserved: SMLib's third span is
empty, so every point the solver finds (all carry tree nodes) skips it.

A curve run along the unit plane's edge is listed as `grefuse`: the
rebuild's boundary search meets the plane's straight edges, and SMLib's
point search on a line (`SmLine`) is analytic, which the port refuses.
SMLib's answer is in the raw evidence.

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
