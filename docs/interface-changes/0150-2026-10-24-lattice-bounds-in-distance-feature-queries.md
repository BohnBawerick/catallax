## 2026-10-24: lattice bounds in distance-feature queries

Pieces: field and query-only comparison
Change: Periodic lattice leaves report no finite bound to the distance-feature candidate visitor, matching the general candidate and display-bound policies. CSG propagation recovers finite bounds from the clipping region.
Why: The lattice definition is unbounded before clipping; its empty analytic constraint list is not a bounded convex cell.
Migration: No public signature, part-construction, point-classification, geometry-limit or exchange-policy change. Candidate generation remains an explicit query and retains caller cancellation.

