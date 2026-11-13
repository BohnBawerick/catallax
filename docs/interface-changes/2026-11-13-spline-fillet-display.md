# Certified display of a spline fillet body

Pieces: shapes, mesh

`mesh_body` now prepares a checked cut chart for one spline fillet between two
polynomial bilinear patch faces with planar ends. Its trimmed side faces use
exact ribbon parameterizations. The existing shared-station mesher also covers
the remaining complete patch faces, including single-chart atlases.

The fillet grid follows the two cap-cut graphs. Complete derivative boxes bound
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
