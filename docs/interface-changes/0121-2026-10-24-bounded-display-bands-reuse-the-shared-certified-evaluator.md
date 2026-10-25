## 2026-10-24: bounded display bands reuse the shared certified evaluator

Pieces: display mesh

Change: a bounded torus band restricts its face-owned chart in place and uses the shared outward-rounded display evaluator with the existing ball fallback. This supersedes the separate band-local frame and angular cache described in the preceding exact-fillet display entry. The merged shared-rail ear guard and cross-face interior-chord ownership check are both retained.

Why: integrating the curved-display work must not bypass its certified sampler or duplicate cached definitions and chart copies.

Migration: no public API, exact geometry, face or edge IDs, default 0.01 mm tolerance, vertex limit, arithmetic allowance, manifest or STEP changes.

