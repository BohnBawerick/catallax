## 2026-11-12: spline fillet faces as body surfaces

Pieces: 2
Change: `SplineSupport` gains `SplineFillet`, a fitted rolling-ball face between two spline faces. `SplineFillet` gains `transformed`, `shares_definition_nodes`, `supports`, `rails`, `enclosed_controls`, `parameters_of_box`, `on_boundary`, `normal_at`, `side_at` and `line_hits`.
Why: the face must stand on a body, like the other fitted rolling-ball faces ([prerequisite](2026-11-12-spline-fillet-caps.md)).
Migration: none. `SplineSupport` is non-exhaustive, so every match outside the shapes crate already has a wildcard, and for the new face those sites keep their refusals.
