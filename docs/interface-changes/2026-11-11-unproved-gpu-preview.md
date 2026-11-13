# Unproved GPU preview

F9 adds `Plan::unproved_preview_wgsl() -> Result<String, Refusal>`. It translates supported compiled field leaves and Boolean instructions to a WGSL `field(p: vec3<f32>) -> f32` function. The result is an unproved display program with no error bound. Its values must never feed an export, a proven measurement or a proof. Existing deterministic evaluation methods and their work accounting are unchanged.

The translator supports stock solids, finite and through extrusions, full and partial revolutions, analytic cells with plane, cylinder, cone, sphere and torus supports, and admitted planar, polyhedral, round, annular and oblique draft cells. It uses the existing processor caches for geometry, angular sectors and through-cut membership. Unsupported leaves and constants outside binary32 range retain typed refusals. An unsupported translation names the leaf index, total leaf count and that leaf's result-face count. A refusal returns no partial program. The existing hardware-preview fence permits the new source generator and the standalone viewer only.

The `catallax-preview` binary consumes this program through wgpu. A missing GPU adapter or unavailable device uses the processor preview. The binary renders P3 and a junction-radius edit, checks the outlines against independent B-rep masks, and writes an unproved named-face ID per pixel on the GPU. A P3 specialization queries the thirteen closed analytic patches; 1,024 seeded surface pixels per frame are checked independently through exact F8 after readback. Both answers and the pixel are retained for every mismatch or reference refusal. This prototype has no SolidDraft integration and no export route. Its [usage and artifact formats](../../crates/catallax-preview/README.md) describe the review output.


## General named preview picking

`catallax_preview::picking::Picker` takes a `catallax_field::display::PreviewMesh`
and snapshots its named triangles. `Picker::points` returns tentative triangle
candidates and triangle-distance estimates in input-point order. Both backends
round the picture and query coordinates to binary32; the processor then measures
in binary64 and the graphics card in binary32. Each point visits triangles in
input order. Equal minima select the first triangle. These are unproved
candidates, not accepted face names.

`Picks::resolve(prepared, points, options, dec)` resolves each input through exact
F8 at its binary32 coordinates. `ResolvedPicks::face(index)` exposes only the
proved owner. This corrects candidates near approximated trim borders. An exact
ambiguity has no accepted handle and retains `Ownership::Refused` with its
`field.ownership-ambiguous` code, contender bounds, spread and cell count.
Proximity to a selected face never licenses accepting its name. Coordinate-count
mismatches refuse with `preview.ownership-points` and both counts. Exhausted
resolution keeps the original numeric refusal and its undecided evidence, with
the point index, total count and coordinates. All queries share the caller's
active allowance and decider; resolution never starts or enlarges an allowance.
The caller supplies `Prepared` from the same revision as the picture.

The picture still carries no geometric proof. Missing chart cells and triangles
that collapse in binary32 stay visible as counts on every batch. The picker uses
no B-rep mesh and no fixed P3 dimensions. It needs the named chart preview; the
coarse surface-nets first paint remains unnamed. F8 supplies all accepted
ownership and its border-margin bound. No preview distance establishes an owner.

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


The planar draft translation omits the projected foot's own plane from its
trim test. The foot lies on that plane by construction. Testing it again can
round it outside the plane and select a farther edge instead. The processor's
unproved draft evaluator uses the same rule. Two P1 draft controls have exact
expected distances `y - 5` to the flange wall; they exercise the binary32 and
binary64 failures separately. This changes tentative draft magnitudes only.
The exact evaluator and its work accounting are unchanged.


## Qualification of 13 November 2026

The required-device run used an NVIDIA GeForce RTX 2070 Max-Q through Vulkan.
All 17 shader cases evaluated the same 20,000 seeded binary32 points. There were
zero sign differences from the processor over 340,000 points, with no
boundary-distance exclusion. The largest value difference was 0.000009850 mm. Repeated GPU
runs and dispatches split into 1,003 points returned identical bits. Processor
values matched bit for bit at 1, 2 and 4 workers. Every generated program also
passed Naga validation without a device or proof work.

The named chart previews used a requested 0.5 mm visual setting under the
unchanged 20M application allowance. That setting is not an ownership tolerance.
Sixteen triangle-centroid samples per part went through the processor and GPU
candidate routes, then `Picks::resolve` and an independent exact F8 query at the
same binary32 coordinates. Every proved reference sample required exact owner
name equality. Every ambiguous reference sample required the identical typed
refusal record and no accepted handle. No proximity test admits a name.

| Part | Triangles | Chart work | Exact accepted names / proved references | Typed refusals | Raw candidates corrected | Accepted mismatches | Resolution work |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| P1 | 6,098 | 781,024 | 14 / 14 | 2 | 0 | 0 | 1,318,736 |
| P2 | 1,152 | 70,316 | 16 / 16 | 0 | 0 | 0 | 714,274 |
| P3 | 2,162 | 54,481 | 15 / 15 | 1 | 1 | 0 | 2,246,201 |
| P4 | 4,360 | 36,285 | 16 / 16 | 0 | 1 | 0 | 511,280 |

All 61 proved reference owners matched. The three refusals retain their
`field.ownership-ambiguous` contender distance bounds, spread and cells. Each
part's 16-point resolution shared one existing 100M verification allowance,
without resets or increases. The two raw candidate disagreements near curved
trims were corrected from F8 proof before exposing accepted handles.

The four pictures had zero missing cells and zero collapsed triangles. GPU and
processor candidates agreed at all 64 samples. Processor candidates and distances
matched at 1, 2 and 4 workers, as did repeated GPU candidates. Separate exact
resolution controls returned identical owner/refusal records and charged work at
1, 2 and 4 workers. Empty resolution charged zero work under a zero allowance;
a 1,000-unit allowance retained `numbers.operation-budget` with the stopping
point and total count. The existing F9 P3 closed-patch route is unchanged.
The picture remains unproved even when F8 proves an ownership query.

Authored P5 radius 2 still refuses compilation with
`fillet.unsupported-overflow`, after 3,367 work in this qualification. No shader
or picker can claim that refused revision. The radius-1/2 edit passes shader
qualification as a separate recipe. This extension does not qualify SolidDraft
first-display integration or M26 timings. Raw evidence retained privately.
