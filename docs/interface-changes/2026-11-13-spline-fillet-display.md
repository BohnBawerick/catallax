# Certified display of a spline fillet body

Pieces: shapes, mesh

`mesh_body` now prepares a checked cut chart for one spline fillet between two
polynomial bilinear patch faces with planar ends. Its trimmed side faces use
exact ribbon parameterizations. The existing shared-station mesher also covers
the remaining complete patch faces, including single-chart atlases.

The fillet grid follows the two cap-cut graphs. Derivative queries reuse exact
Bézier span views translated by a constant origin; translation preserves every
derivative and avoids interval cancellation against world coordinates. Complete derivative boxes bound
the barycentric chord error across carrier knots. A fixed projection proves
positive symmetric secant matrices for every cell and includes output rounding.
The cut graph's parameter error contributes to the surface bound. Shared cut
vertices lie exactly on their cap planes, and their displacement contributes to
the chord and embedding checks. The input body and its fit marks remain intact.

`CutCurve.parameter_upper` reports the uniform error in the linear graph's tau
coordinate before spatial-control rounding. The STEP curve, its spatial bound,
and its output controls remain unchanged. Callers constructing `CutCurve`
literals must supply this bound.

Unsupported trims or unproved grids return `mesh.spline-fillet-unproven` with
the requested tolerance and vertex cap, plus any available failing bound.
Existing precision and resource refusals remain typed. No mesh is returned when
the 0.01 mm default bound, shared boundary, embedding, or closure check fails.

`RationalPatch::display_grid_adaptive` is an opt-in grid entry that retains
prepared span and cell bounds, caches accepted cells, and bisects only a failed
cell's limiting interval. It preserves exact knots and shared stations. The
existing `display_grid` entry retains its refinement and output contract.

`BSplineSurface::derivatives_box` encloses first derivatives over an exact
rational box. It evaluates the prepared derivative nets on each intersected
owned knot span, proves each positive quotient, and joins the enclosures.
At a knot it includes both one-sided derivatives. Existing binary-box entries
keep their behavior.
