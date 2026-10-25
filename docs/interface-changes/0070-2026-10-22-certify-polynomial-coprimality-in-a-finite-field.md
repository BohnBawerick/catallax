## 2026-10-22: certify polynomial coprimality in a finite field

Piece: 1
Change: polynomial GCD first tries a fixed 31-bit prime. It concludes coprimality over the rationals only when every coefficient denominator is invertible, both degrees are preserved, and finite-field Euclid ends at a nonzero constant. Otherwise the existing exact polynomial GCD runs unchanged.
Why: square-free annihilators repeatedly entered large integer remainder and content calculations merely to establish that the GCD with their derivative was one.
Migration: no public signatures, coefficient domains, default budgets or deadline changes. Normalized GCDs stay unchanged. Regressions cover a common factor lost by modular degree reduction, noninvertible rational denominators and coprime rational polynomials whose modular images share a factor.

