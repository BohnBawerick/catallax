## 2026-10-22: sampled two-kernel display comparison

Pieces: 5, field queries, recipe consumers, display comparison
Change: add `catallax-compare` with `Controller`, revision tokens, scoped `Options`, complete `Report` values and explicit agreement/disagreement/unresolved outcomes. Add piece 5 query helpers `face_contains` and `surface_witness`, and `Plan::is_proven_empty` in the independent field evaluator. No existing geometry builder, recipe schema, question, manifest or STEP signature changes.
Why: the background display check needs certified distances in both directions, exact trimmed-face associations, face origin sets and cancellation that cannot publish stale or partial work.
Migration: existing callers are unchanged. Comparison callers begin a revision per edit and dispatch its check on their own worker. Read `Controller::latest` only as a revision-tagged snapshot. Treat regional and sampled agreement as scoped display evidence, never as a 1 nm or whole-boundary certificate. `surface_witness` supplies an exact rational chart projection, not a closest point; measure its distance and call `face_contains` before using it as an upper-distance witness.

