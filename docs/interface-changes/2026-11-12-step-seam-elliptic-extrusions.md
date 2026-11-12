## 2026-11-12: oblique extrusions meet ND12

Pieces: 7, judge interface, verification exchange

Change: under `Seams::Required` (`exact-brep-question/3`, `catallax-step/1`) the elliptic wall of an oblique extrusion (`SURFACE_OF_LINEAR_EXTRUSION` over an `ELLIPSE`) gets a chart and a pcurve for each edge, where [the seam contract](2026-11-11-step-seam-contract.md) refused it with `step.seam-unsupported`. The wall is charted like a cylinder's: `u` is the swept ellipse's own parameter and `v` is the height along the extrusion axis. The seam is the ruling at the ellipse's major vertex or at its opposite, the only places an ellipse can start its parameter. A wall that crosses the seam is cut there and the seam is written as a `SEAM_CURVE`.

A wall whose boundary lines lie along the seam is charted without a cut. The seam is then its own boundary: a half wall of a slot has its two end rulings on the major axis, so one ruling is the seam at `u = 0` and the other is interior, or the first sits at `u = 2 pi`. The chart gives the vertices of that ruling the parameter of the side the face lies on, and the ruling's pcurve is fixed at it. A wall that would need both a cut and a boundary on the seam keeps its refusal.

Why: PLAN M24 asks for seam curves in STEP. Two of the 19 public recipes that refused under ND12 are oblique extrusions of a circle and of a slot.

Measured, raw evidence retained privately: `algebraic-extrusions/circle` and `algebraic-extrusions/slot` write under ND12. Each file reads in the pinned `cadquery-ocp==8.0.1.1.0` with no ShapeFix, SameParameter or tolerance limiting, passes `BRepCheck_Analyzer`, and has the certified volume to better than 1e-9 relative. The ND12 check agrees. Each file has the same bytes at 1, 2 and 4 kernel threads. The third recipe of that fixture directory, `rectangle`, had no elliptic wall and is unchanged. Seven of the 195 public STEP files still refuse under ND12; all are offset parts.

Version 2 output is unchanged: the new routes are reached only under `Seams::Required`, and the 227 cases of the 226 public recipe files written without the wrapper have the manifests and bytes they had before. Under `Seams::Required` only the two files above change. No saved replay answer, tolerance or work allowance changes.

Migration: none. Callers that handled `step.seam-unsupported` for these shapes now receive a file; callers of `exact-brep-question/2` see no change.
