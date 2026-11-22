## 2026-11-22: field spline sweep ranges

Pieces: field, mesh and compare
Change: `Spine::prove_spline_plane` accepts the offset of a plane chart.
`spline_spine::sweep::Sweep` bounds the distance to the finite short-arc
face, with optional proved cap cuts. `Sweep::distance_at_most` refines
the range against a distance bound. The rounded shell classifies solid
membership. A rounded loft exposes its blend origin and guaranteed
blend range through `Plan::spline_blend_origin` and
`Plan::spline_blend_distance`, which refines against the supplied bound.
`catallax_mesh::display_spline_fillet_face` supplies a measured B-rep
display of the selected trimmed face at the requested chord tolerance.
`display_spline_fillet_carrier` supplies the positional grid and bound
for an entire fitted carrier without body cuts.
Why: a field fillet should prove its own sweep and compare it against the
B-rep without depending on the B-rep's offset or fitted carrier.
Migration: callers of the spline sweep set the part tolerance when proving
its spine. See [spline fillets in comparison](../FIELD.md#spline-fillets-in-comparison)
for supported recipe inputs and copy behavior.
See [selected spline fillet comparison](../COMPARE.md#selected-spline-fillet)
for report semantics. Sweep refusals are reported under
`field.spline-sweep`, `field.spline-sweep-reach`,
`field.spline-sweep-sector`, `field.spline-sweep-cap`, or
`field.spline-sweep-clearance` as applicable.
