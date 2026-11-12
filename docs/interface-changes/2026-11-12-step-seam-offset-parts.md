## 2026-11-12: offset parts meet ND12

Pieces: 7, judge interface, verification exchange

Change: under `Seams::Required` (`exact-brep-question/3`, `catallax-step/1`) offset and distance-limited parts that mix cone, cylinder, sphere and torus faces get charts, seams and pcurves on every curved face, where [the seam contract](2026-11-11-step-seam-contract.md) refused them with `step.seam-unsupported`. Four export-only steps on a copy of the model make the existing chart writers reach them:

- A torus band that ends on the outer equator from below, such as the lower round of an offset cylinder, is written whole on the same torus with its axis reversed. Its tube angle then starts on the equator like the upper band's, which is the only start the chart accepts. Bands across the equator keep the cut of [trimmed torus charts](2026-10-27-trimmed-torus-step-charts.md).
- A spherical face whose pole the model marks with a vertex loop loses that loop, and its pole cap, which the chart writers already add for an interior pole, carries the pole. The written file has no `VERTEX_LOOP` there; ND12 still admits one on a sphere.
- A vertexless torus ring that the model stores as an intersection component but that is a circle, as a coaxial cylinder or plane section is, is replaced by that circle, with the same incidences. The ring's `CIRCLE` entity and its bytes are unchanged; only the chart's cut queries stop treating it as a quartic.
- A point on the written outer equator that rounds a hair below it takes the tube angle `-epsilon`, not a turn above, when it is within one millionth of a turn of the seam. A fit that mixes the two lifts has an error of a whole turn and is refused by the deviation check.

Why: PLAN M24 asks for seam curves in STEP. Seven of the 19 public recipes that refused under ND12 were offset or mixed parts.

Measured, raw evidence retained privately: `cylinder-distance/coaxial-cylinders`, `sphere-cylinder`, `sphere-back` and `sphere-front`, `cone-distance/tip-outward-posed`, and `curved-limits/torus-entry-next` now write under ND12. Each file reads in the pinned `cadquery-ocp==8.0.1.1.0` with no ShapeFix, SameParameter or tolerance limiting, passes `BRepCheck_Analyzer`, and has the certified volume to better than 1e-9 relative (all six within 1e-13). The ND12 check agrees. Each file has the same bytes at 1, 2 and 4 kernel threads.

One recipe stays refused, with its reason: `cone-distance/tip-outward-small`. Its offset sphere has radius 1/200 at the tip of a cone 8 long, so a trim ring lies so near the cone's axis that the chart's fit needs more spans than its work allows. It refuses `step.budget-exhausted` ("UV chart domain cannot be separated from its axis"), as it has since the cone charts reached it, and the allowance is not raised. Of the 195 public STEP files, one refuses under ND12.

Version 2 output is unchanged: the new routes are reached only under `Seams::Required`, and the 227 cases of the 226 public recipe files written without the wrapper have the manifests and bytes they had before. Under `Seams::Required` only the six files above change from the previous interface note. No saved replay answer, tolerance or work allowance changes.

Migration: none. Callers that handled `step.seam-unsupported` for these shapes now receive a file; callers of `exact-brep-question/2` see no change.
