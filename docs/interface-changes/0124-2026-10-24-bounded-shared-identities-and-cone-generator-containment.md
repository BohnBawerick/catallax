## 2026-10-24: bounded shared identities and cone generator containment

Pieces: 1, 4, 5
Change: `Decider::rational_identity(&Real)` optionally proves a rational constant using only bounded rational-function identities over shared definitions, after domain validation. An unrecorded complete line on a cone generator now refuses with `intersect.partial-containment-required`, because only half the line belongs to the retained nappe.
Why: circular-section adjacency needs cheap exact zero proofs before interval refinement, while full generator containment would discard the cone domain.
Migration: callers may use the optional identity before existing certified fallbacks; `None` proves nothing. Curve/surface consumers must preserve the new named refusal when partial containment is required. Existing recorded incidence remains authoritative. No budget, tolerance or file format changes.

