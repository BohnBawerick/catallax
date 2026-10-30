## 2026-10-31: finite straight edge-pair readout

Pieces: 6 (measurement), SolidDraft Catallax adapter

Change: Added `catallax_measure::edge_pair` and `EdgePairRequest`. For two finite line edges with exact rational ends, it returns the minimum distance over both closed segments, the distance between their native half-arc-length midpoints, or the unoriented direction angle in degrees. Rational distances are exact; irrational roots are outward enclosures of nonzero width. The angle answers only exact parallel (0) or perpendicular (90) carriers. New refusal codes are `measure.irrational-endpoint`, `measure.degenerate-edge` and `measure.angle-domain`; curved or unbounded edges reuse `measure.unsupported-curve`.

Why: SolidDraft V2 batch 2 measures finite edge pairs. The crate had edge length and midpoint, face and point distances, but no certified edge/edge readout.

Migration: Existing callers need no changes. The adapter resolves persistent references to current edge handles first, keeps its own handle errors for unknown or stale references, and maps the native refusal codes as causes. Work exhaustion keeps its native code (`numbers.operation-budget` or `numbers.budget-exhausted`).
