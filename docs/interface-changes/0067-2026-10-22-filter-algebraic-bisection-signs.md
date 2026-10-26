## 2026-10-22: filter algebraic bisection signs

Piece: 1
Change: when the odd part of a bisection denominator exceeds 512 bits, root enclosure tries a budgeted polynomial ball at the caller's working precision, capped at 512 bits, before exact rational Horner evaluation. Smaller odd denominators retain the direct exact path. Only a strict interval sign selects a side; zero-containing enclosures retain exact evaluation. Bisection points, root selection and final rounding stay unchanged.
Why: a certified bracket can have very large odd denominators without needing equally large working precision. Exact polynomial values at every midpoint exhausted the budget even for the real root of `x^3 - 2`.
Migration: no public signature, default allowance or deadline changes. A regression brackets `cbrt(2)` with 4,096-bit odd denominators and requires width `1e-30` within one million work units; GMP cubic inequalities check containment. The cold/warm precision-history matrix checks work identity and insufficient-allowance fallback.

