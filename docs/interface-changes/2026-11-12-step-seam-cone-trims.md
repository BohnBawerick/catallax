## 2026-11-12: trimmed cones and open conics meet ND12

Pieces: 7, judge interface, verification exchange

Change: under `Seams::Required` (`exact-brep-question/3`, `catallax-step/1`) a cone face with any closed trims now gets a chart, a seam and a pcurve for each edge, where [the seam contract](2026-11-11-step-seam-contract.md) refused it with `step.seam-unsupported`. The chart is the cylinder's: the face is cut at a meridian, each edge carries a certified degree-seven pcurve, and the seam is a `SEAM_CURVE`. Plane sections of the nappe that are parabolas or hyperbolas now write a pcurve too; the 3D `PARABOLA` and `HYPERBOLA` entities are unchanged, and the pcurve shares their parameter.

A cone face is charted only when the apex is proved outside it: no trim passes through the apex, and the crossings of the chart meridian alternate from entering to leaving in height order. A cone face with the apex on a trim or inside the face keeps its refusal. A scene with any other face outside the chart writers also keeps it.

Why: PLAN M24 asks for seam curves in STEP. The 19 public recipes that refused under ND12 included eight cones trimmed by planes or limits, and two offset parts whose only unseamed face was a cone.

Measured, raw evidence retained privately: the eight cone recipes (`open-conics/hyperbola` and `parabola`, `cone-distance/oblique-inward`, `curved-limits/cone-entry-next` and `cone-entry-to`, and the three `polygon-limits/cone-*` files) and the two offset recipes `cone-distance/tip-outward` and `tip-outward-rational` now write under ND12. Each file reads in the pinned `cadquery-ocp==8.0.1.1.0` with no ShapeFix, SameParameter or tolerance limiting, passes `BRepCheck_Analyzer`, and has the certified volume to better than 1e-9 relative. The ND12 check agrees on all ten. Each file has the same bytes at 1, 2 and 4 kernel threads. Nine of the 195 public STEP files still refuse under ND12; they are the two elliptic extrusions and seven offset parts.

Version 2 output is unchanged: the new routes are reached only under `Seams::Required`, and the 227 cases of the 226 public recipe files written without the wrapper have the manifests and bytes they had before. Under `Seams::Required` the other 216 cases give the same manifests and bytes as before, apart from the ten files above and `cone-distance/tip-outward-small`, whose refusal moves from `step.seam-unsupported` to `step.budget-exhausted` because its chart attempt now runs and exhausts its precision. No saved replay answer, tolerance or work allowance changes.

Migration: none. Callers that handled `step.seam-unsupported` for these shapes now receive a file; callers of `exact-brep-question/2` see no change.
