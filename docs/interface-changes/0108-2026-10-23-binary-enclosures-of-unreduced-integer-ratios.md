## 2026-10-23: binary enclosures of unreduced integer ratios

Pieces: 1, 5
Change: `Ball::from_ratio(&Integer, &Integer, u32) -> Option<Ball>` encloses an integer ratio without constructing a normalized `Rational`. Either denominator sign is accepted; a zero denominator returns `None`.
Why: cylindrical winding repeatedly normalized exact diamond-chart corner fractions before immediately rounding them to binary balls. Homogeneous corner comparisons and direct ratio enclosures avoid those GCDs while preserving the exact enclosure.
Migration: existing callers are unchanged. Piece 5 uses the new constructor for phase extrema. No budget, tolerance, answer format or refusal code changes.

