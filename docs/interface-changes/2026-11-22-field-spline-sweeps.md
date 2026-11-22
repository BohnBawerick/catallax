## 2026-11-22: field spline sweep ranges

Pieces: field, mesh and compare
Change: `Spine::prove_spline_plane` accepts the offset of a plane chart.
`spline_spine::sweep::Sweep` bounds the distance to the finite short-arc
face, with optional proved cap cuts. Its `query` also classifies membership
of the swept balls. A rounded loft exposes its blend origin and guaranteed
blend range through `Plan::spline_blend_origin` and
`Plan::spline_blend_distance`. `Plan::rounded_loft_face_bound` supplies an
upper distance from a certified retained side or cap point.
`catallax_mesh::compare_position_mesh` supplies a measured B-rep display
with a complete positional bound at the requested chord tolerance.
Why: a field fillet should prove its own sweep and compare it against the
B-rep without depending on the B-rep's offset or fitted carrier.
Migration: callers of the spline sweep set the part tolerance when proving
its spine. A rounded loft builds through the ordinary recipe compiler.
Compare uses the positional display for the named spline loft blend and
reports missing reverse coverage. The positional display carries no
facet normal or export certificate. Refusals are reported under
`field.spline-sweep`, `field.spline-sweep-reach`,
`field.spline-sweep-sector`, `field.spline-sweep-cap`, or
`field.spline-sweep-clearance` as applicable.
