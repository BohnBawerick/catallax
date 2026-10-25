## 2026-10-22: retain owned multivariate accumulators

Piece: 1
Change: owned polynomial addition retains its left term tree and moves the right terms into it. Entry-based accumulation moves vacant coefficients directly, computes occupied sums once, and removes exact cancellations.
Why: curve substitution repeatedly cloned growing term trees and copied incoming arbitrary-precision coefficients before discarding them.
Migration: operator signatures, coefficient domains, term order and normalized results are unchanged. Borrowed addition retains its independent result. Targeted independent polynomial tests cover substitution, cancellation, derivatives, resultants and certified systems; no precision or work-policy changes.

