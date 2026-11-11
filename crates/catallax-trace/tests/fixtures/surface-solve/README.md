# Surface point solver gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for `SmSurface::LocalPointSolve_1`,
the local surface point solver that `SmSurface::LocalPointSolve` runs for
minimising, normalising and intersecting (the found flag, the parameters, the
distance, the termination reason and the iteration count), and for
`SmSurface::GlobalPointSolve_1`, the global solver behind
`SmSurface::GlobalPointSolve` (every solution's parameters and distance). The library build is
the one described in [`../foundation/README.md`](../foundation/README.md).
Catallax never compiles or links this C++ code.

The surfaces are those of the surface evaluation gate (a bicubic bump with
interior knots, polynomial and rational, a rational loft-like patch, a cubic
patch whose u = 0 edge collapses to a pole, and a bilinear patch) and a
rational cylinder closed in u. The cases drop points from above and below,
from beyond every side and corner (the answer stops on the boundary through
SMLib's iso-curve boundary solution), onto subdomains, from points on the
surface, from interior knots and corners, across the cylinder's seam, from
the cylinder's axis, where every direction is equally near, from deep inside
the cylinder, where SMLib limits the second-derivative terms of the Newton
step, and from the pole, where SMLib inverts the point near the singularity
and finishes by dropping it on the iso-curve through the pole from several
guesses, minimising, normalising and intersecting.
A fixed pseudo-random sweep of 80 test points and guesses covers the five
surfaces.

The `global` cases solve from points above, beside and beyond the surfaces, on
them, at the pole, on the cylinder's axis and just across its seam (where the
answer is copied to the other side of the seam), on a rational fan whose first
row of control points lies on a circle about the test point (every point of
that row ties for nearest), and from a sweep of 40 more points. They reach the
control-net triangle guess, the control-point guess with its leg refinement
and ties, the shrinking parameter grids with their stride check, the plane-step
refinement, the second grid search after a failed local solve, and the seam
copies. The surfaces are built as SMLib's `CreateCanonical` builds them: a net
whose weights are all equal is stored as polynomial.

Binary64 values match within 16 units in the last place; every one matched bit
for bit when the fixture was made. Mutating the step limit, the line-search
fractions, the bracketed search tolerance, the derivative floor or the
second-derivative limit of the local solver, or the grid size, plane-step limit
or stride limit of the global solver, changes cases. The Newton improvement
ratio and the boundary-hit counts of the local solver only act in windows these
cases do not reach, and several global heuristics (the net-shape limits, the
near-net-point test, the seam tolerance, the number of grid passes, the
curvature noise limit, the early stop among ties) change the guess but not the
local solver's answer from it in these cases.

To regenerate, build the driver as in the foundation README with
`crates/catallax-trace/tests/fixtures/surface-solve/oracle.cpp`, then run

```sh
python3 crates/catallax-trace/tests/fixtures/surface-solve/generate.py "$ORACLE" "$RAW_DIR"
```

Raw evidence retained privately.
