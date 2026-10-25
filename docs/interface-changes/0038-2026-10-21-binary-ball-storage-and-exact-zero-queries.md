## 2026-10-21: binary ball storage and exact-zero queries

Pieces: 1, 4, 6
Change: piece 1 adds `Ball::is_exact_zero(&self) -> bool`, true only for the singleton `[0, 0]`, without allocating field copies. Ball midpoints and radii use private odd integer mantissas with binary exponents. Reciprocal rounding uses integer quotients of dyadics instead of normalizing rational endpoint fractions. Rational accessors, rounding, enclosure formulas and work charges are unchanged.
Why: validated Taylor arithmetic spent most of its time allocating and normalizing rational binary denominators, including copies made only to test a coefficient for zero. Binary exponents and exact quotient rounding remove those denominators and their repeated GCDs.
Migration: trace and measurement coefficient checks use `is_exact_zero`; the allocating helper functions are removed. Public rational accessors retain their signatures and canonical values. A before/after arithmetic probe retains identical midpoint/radius hashes, and the GMP/MPFR enclosure suites check containment independently.

