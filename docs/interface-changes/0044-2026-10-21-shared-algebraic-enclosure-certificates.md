## 2026-10-21: shared algebraic enclosure certificates

Piece: 1
Change: immutable irrational-root snapshots share their polynomial, interval and a bounded memo of eight working-precision enclosures. A cache hit returns the exact ball produced by the original bisection at that precision and replays its work cost. An insufficient allowance takes ordinary computation to preserve partial-work accounting. Refinement and defining-polynomial replacement do not inherit another snapshot's memo.
Why: adjacency coefficient checks repeatedly enclosed the same algebraic leaves. Public root-list memoization did not retain those bisection certificates. Snapshot sharing also avoids copying large polynomial coefficients and interval endpoints whenever an algebraic value is cloned.
Migration: no public signature, default limit, rounding rule or algebraic-equality contract changes. No higher-precision substitution occurs. Independent GMP checks cover budget-limited cold/warm results and a cached request still refuses an expired deadline.

