## 2026-11-12: a fillet between two spline faces on a body

Pieces: shapes, topology, features
Change: a `fillet` step that selects one convex edge between two one-span patch faces now rounds it with a fitted `SplineFillet` face, trims both faces along the rails and cuts each end plane along the face. It adds the refusal code `fillet.unsupported-spline` for a concave edge, a side face that is not its whole patch, or a fitted face that would face away from the ball. `SplineFillet` gains `rail_parameters`, `rail_point`, `rail_within`, `rail_curve` and `faces_ball`. `RationalPatchFace` gains `of_chart_face`, which also checks a face on an atlas of one chart as that chart's patch.
Why: M11 needs the fillet face on a body, with stable names under a radius edit ([prerequisite](2026-11-12-spline-fillet-surfaces.md)).
Migration: such a selection used to refuse with `fillet.unsupported-contact`. It now builds, or refuses with a fillet or spline-fillet code naming what was not proved.
