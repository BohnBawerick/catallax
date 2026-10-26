## 2026-10-22: try separation certificates before contact proofs

Piece: 5
Change: adjacency probes try boundary equations with zero-excluding cached coefficients first. Endpoint checks likewise try coordinates whose source and endpoint enclosures are disjoint before potentially coincident coordinates.
Why: an exact contact in the first coordinate could trigger expensive algebraic elimination even when another coordinate immediately separated the whole path.
Migration: cached bounds select only the order of existential checks. The same exact whole-path and selected-component predicates still decide adjacency; no public signature, default allowance or deadline changes.

