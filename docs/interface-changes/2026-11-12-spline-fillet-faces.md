## 2026-11-12: rolling-ball fillet faces between spline faces

Pieces: 2
Change: `catallax-shapes` adds `SplineFilletDefinition` and `SplineFillet`, a checked rolling-ball face between two single-span rational patches, with its rails as `RailSpan` lists and the target constant `SPLINE_FILLET_TARGET_NM`. It adds six refusal codes: `fillet.spline-definition`, `fillet.spline-deviation`, `fillet.spline-spine`, `fillet.spline-support-spans`, `fillet.spline-clearance` and `fillet.spline-carrier-irregular`.
Why: M11 needs a proved fillet face on an edge between two spline faces, on the spine that the exact offset meeting proves ([prerequisite](2026-11-12-exact-offset-patch-meetings.md)).
Migration: none; the change is additive and nothing calls it yet.
