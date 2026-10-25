## 2026-10-20: real roots of degree up to 4 carry irreducible defining polynomials

Pieces: 1 (and every piece that calls `real_roots`)
Change: a root `real_roots` returns that is not rational carries a defining polynomial irreducible over Q whenever its square-free factor has degree at most 4, so a quartic that splits into two rational quadratics is split. Piece 1's case 40 gains such a quartic and an exactly zero expression on its roots that must be `ProvenCoincident` under the default budget. No signature, outcome, code or budget default changes.
Why: piece 4's hits of the case 29 bite's branch with the plane `z = 0` are roots of `(x^2 - a)(x^2 - b)` with rational `a`, `b`. Kept as roots of the quartic, an exactly zero square-root argument at those hits could not be proven within 10^9 work units, so `intersect_curve_surface` refused `numbers.sqrt-not-separated`, and so did piece 5's `classify_point` of a case 55 point with `z = 0`. With the quartic split into its quadratic factors, the same intersection answered with its two crossing hits using 7.8e6 work units.
Migration: piece 1 factors square-free factors of degree 4 with no rational root into rational quadratics where they exist (through the resolvent cubic, for example); rational roots are already removed. Callers need no change.

