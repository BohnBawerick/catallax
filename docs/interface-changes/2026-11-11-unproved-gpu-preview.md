# Unproved GPU preview

F9 adds `Plan::unproved_preview_wgsl() -> Result<String, Refusal>`. It translates supported compiled field leaves and Boolean instructions to a WGSL `field(p: vec3<f32>) -> f32` function. The result is an unproved display program with no error bound. Its values must never feed an export, a proven measurement or a proof. Existing deterministic evaluation methods and their work accounting are unchanged.

The translator supports stock solids, finite and through extrusions, full and partial revolutions, analytic cells with plane, cylinder, cone, sphere and torus supports, and admitted planar, polyhedral, round, annular and oblique draft cells. It uses the existing processor caches for geometry, angular sectors and through-cut membership. Unsupported leaves and constants outside binary32 range retain typed refusals. An unsupported translation names the leaf index, total leaf count and that leaf's result-face count. A refusal returns no partial program. The existing hardware-preview fence permits the new source generator and the standalone viewer only.

The `catallax-preview` binary consumes this program through wgpu. A missing GPU adapter or unavailable device uses the processor preview. The binary renders P3 and a junction-radius edit, checks the outlines against independent B-rep masks, and writes an unproved named-face ID per pixel on the GPU. A P3 specialization queries the thirteen closed analytic patches; 1,024 seeded surface pixels per frame are checked independently through exact F8 after readback. Both answers and the pixel are retained for every mismatch or reference refusal. This prototype has no SolidDraft integration and no export route. Its [usage and artifact formats](../../crates/catallax-preview/README.md) describe the review output.


## General named preview picking

`catallax_preview::picking::Picker` takes a `catallax_field::display::PreviewMesh`
and snapshots its named triangles. `Picker::points` returns one unproved face
handle and triangle-distance estimate per input point. The caller resolves each
handle against the same revision's field registry. Both backends round the
picture and query coordinates to binary32; the processor then measures in
binary64 and the graphics card in binary32. Each point visits triangles in input
order. Equal minima select the first triangle. The processor's answers are
identical at 1, 2 and 4 workers. Cross-device bit equality is not claimed.

This query needs the named chart preview. The coarse surface-nets first paint
still has no names. The picker uses no B-rep mesh and has no fixed P3 dimensions.
It can pick any named chart preview that the field supplies. Missing chart cells
and triangles that collapse in binary32 stay visible as counts on every batch.
The picker does not fill those gaps or prove that its name is the nearest defined
face. `Prepared::nearest_face` remains the certified route, including its
contender bounds, border-margin lower bound and typed ambiguity refusals.

The caller supplies a point-triangle comparison limit. An oversized product or
more than 4,194,240 points refuses with `preview.ownership-size` before evaluation.
The same code limits the source to 2,796,202 triangles, within the default
128 MiB storage buffer. Non-finite coordinates or magnitudes above 100,000,000 mm
refuse with `preview.ownership-range` and the point. Invalid indices use
`preview.ownership-mesh`. An empty or fully collapsed mesh uses
`preview.ownership-empty`, with triangle, collapsed-triangle and missing-cell
counts. Zero workers use `preview.ownership-workers`; non-finite results use
`preview.ownership-non-finite`. Device errors retain the existing GPU codes.
These display limits do not change a proof allowance or tolerance.
