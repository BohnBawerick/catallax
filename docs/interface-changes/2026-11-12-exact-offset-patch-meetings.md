## 2026-11-12: Exact offset patch meetings

Pieces: 2

Change: `RationalPatch::offset(distance)` returns a new `OffsetPatch`, the exact offset `P + d N/|N|` at a signed rational distance. `OffsetPatch::offset_meeting` returns the existing `PatchPairIntersection` over both supports' parameters, and its `point` gives the offset point when the first side is an offset. `RationalPatch::patch_meeting` is unchanged.

Why: a rolling-ball fillet's spine is where two offsets meet (M11), and an offset is not a rational surface.

Migration: none; the change is additive. An offset meeting uses the existing refusal codes. A cell whose offset is not proved regular with its support's orientation is never accepted, and at the split limit it refuses with `intersect.spline-pair-unresolved`.
