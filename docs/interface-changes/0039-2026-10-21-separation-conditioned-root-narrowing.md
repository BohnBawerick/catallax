## 2026-10-21: separation-conditioned root narrowing

Piece: 1
Change: root isolation uses the exact quadratic discriminant when rational divisor enumeration is incomplete. Irrational closed-form brackets are reserved for branches closer than the magnitude-scaled Newton guard; well-separated branches retain the ordinary isolator's brackets. Rejected magnitude-based Newton proposals can use interval Newton. Cross-factor separation skips only refinement targets that leave every interval unchanged, retaining the original power-of-two target grid. Exact signs still certify every bracket; the guard selects an exact algorithm, not geometric coincidence.
Why: translated clustered roots defeated magnitude-only proposals. The quadratic smoke changes from default-budget refusal at 999,914,807 work units to two proven irrational roots with 900,624; the analogous cubic changes from refusal at 999,941,208 to one proven irrational root with 214,834,275.
Migration: no public signatures, default budgets, rationality-width targets or irreducibility guarantees change. Independent GMP certificates check the root lists and brackets, including the translated rational cubic.

