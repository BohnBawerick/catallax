## 2026-10-24: ruled display strips share the certified sampling filter

Pieces: display mesh

Change: ruled cylinder and elliptic-cylinder strips certify their existing boundary coordinates through the shared display sampler, including its ball fallback, instead of unconditionally rebuilding ball evaluations. Shared polylines and candidate triangle coordinates are unchanged; conservative reported residual bounds may change.

Why: the new-kind path must reuse the merged certified evaluator rather than introduce avoidable high-precision work.

Migration: no public API, exact geometry, picking IDs, requested controls, tolerance, vertex or work allowances, manifest or STEP changes.

