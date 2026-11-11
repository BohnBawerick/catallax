# Local surface-surface gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for `SmSurface::LocalSurfaceIntersect`
and `SmSurface::LocalSurfaceSolve` (minimising, maximising, normalising): the
found flag and, when found, both parameter pairs and the solution value. The
library build is the one described in
[`../foundation/README.md`](../foundation/README.md). Catallax never compiles or
links this C++ code.

The surfaces are those of the surface point solver gate and canonical planes
over `[-20, 20]²`. The intersection cases meet the bumps with level and tilted
planes from guesses near the meeting curve, with the plane moved each step
(steepest descent) and fixed, with and without side tests and with a previous
point near and far; guesses on the sides of the bumps and of a sub-domain, and
interior guesses whose answer lies beyond a sub-domain side; the loft with
planes; the cylinder with walls along its ruling (tangent), 1e-3, 1e-5 and 1e-7
inside and outside it under tolerances of 1e-5, 1e-4 and 1e-3, and walls where
its normal turns 2 and 4 degrees from theirs; the patch with a pole against
planes through, 1e-9 below and away from the pole; the cylinder with a bump and
a bump with its rational copy; and a seeded sweep of 30 guesses. The solve cases
find nearest, farthest and normal pairs between a bump and planes above and
below it, the cylinder and a bump, the loft and a tilted plane, the two bumps,
and the pole patch and a plane.

Two answers deliberately differ from SMLib. SMLib's side-curve search covers
the whole iso-curve and keeps a meeting, or a final answer, outside the
requested domains; the port keeps both parameter pairs inside them. The driver
records such SMLib answers as no answer and logs them in the raw evidence: four
cases (sub-domains of a bump, and a side curve of a flat patch that meets a wall
only beyond the requested half). And SMLib's Newton solver reports an
evaluator failure at the guess as a found answer, which SMLib returns; the port
returns the typed failure. The unit tests in `tests/surface_surface.rs` cover
both.

Binary64 values match within 16 units in the last place; every one matched bit
for bit when the fixture was made. Changing the 100-to-1 rule for a previous
point, removing the out-of-bounds recovery along the sides, the early answer of
the plane-held evaluator, steepest descent, the plane's origin between the two
points, the side-curve perpendicularity exit, the plane-normal refinement or
its gap shift changes cases. The recursion for an answer that lands on a side,
the iteration and acceptance limits of the plane-held solve, the minimising
refinement under steepest descent, the final normal refinement and the
3-degree near-tangency threshold leave every output unchanged: no case lands
exactly on a side, needs more than 20 iterations or an acceptance between two
and four times the tolerance, and where the refinements run the plane-held or
plane-normal solve has already reached the closer pair.

To regenerate, build the driver as in the foundation README with
`crates/catallax-trace/tests/fixtures/surface-surface/oracle.cpp`, then run

```sh
python3 crates/catallax-trace/tests/fixtures/surface-surface/generate.py "$ORACLE" "$RAW_DIR"
```

Raw evidence retained privately.
