## 2026-10-24: planar draft caps on cylindrical neighbours

Pieces: features, field

Change: `draft::capped::CappedCylinder` retains an unchanged cylinder and two opposing planar cap definitions. Exact affine heights define its disk footprint, strict material, neutral-section area and finite enclosure. A cap intersection can remove part of that footprint without introducing a cap or origin.

`cap_curve` and `FootprintBoundary` construct native circles, ellipses, chord endpoints and tangent contacts. The B-rep assembles their shared topology directly; the provisional enclosing-solid Boolean construction was removed. Tangent cap contact uses one side wire, preserving a manifold convex boundary.

The field publishes an independent bounded analytic cell with residual distance semantics. Native charts survive rigid copies and supply current face candidates by support polynomial. Its strict witness also supports later drafts whose footprint excludes the cylinder axis.

Integration: draft region construction, rigid copies, prepared planar trims and native cap cells propagate the merged fallible analytic-constraint API. Evaluated definitions retain upstream certification snapshots alongside draft consumption history.

Limits: multiple interacting curved cells, other mixed support families and general curved reference disambiguation remain unfinished. The rational candidate interface still refuses genuinely algebraic placements. Successful geometry and STEP import do not imply an interactive comparison-budget pass.

Migration: no new recipe version, face origins, approximation tolerance or work allowance.

