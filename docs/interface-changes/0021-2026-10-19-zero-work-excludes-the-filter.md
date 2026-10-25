## 2026-10-19: zero work excludes the filter

Pieces: 1, 2
Change: clarify that `max_work = 0` permits no decision work, including the f64 filter; a sign returns `NotExact` with a bound from its cached enclosure. Work accounting, defaults, public signatures and implementation behavior are unchanged.
Why: the architect ruled that piece 2's `not_exact_curve_and_point_distance_bounds` test incorrectly expected filter proofs with no work, despite piece 1 charging one work unit for an operation on doubles.
Migration: assertions expecting a filter-proven nonzero coordinate use `max_work = 1`, the smallest sufficient budget; assertions expecting zero-work uncertainty retain `max_work = 0`. Expected outcomes and decision counts stay unchanged.

