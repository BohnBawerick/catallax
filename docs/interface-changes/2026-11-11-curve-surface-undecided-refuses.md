## 2026-11-11: curve-surface intersection refuses an undecided decision

Pieces: 1, 4, 5
Change: `intersect_curve_surface` no longer answers with marks. When a candidate's branch side, a hit's contact or the order of two hits ends `NotExact`, it refuses with that outcome's `numbers.*` code and adds nothing to the store. Hits are always `MeetOf::CurveSurface`; the `MeetOf::Surfaces` form for a hit on an undecided branch is gone. `Mark::from_not_exact` is removed.
Why: each old degenerate branch acted on an undecided answer without a proof. A bound on the distance to a dividing plane is not a distance to the branch, a tangent label changes the crossing parity later stages count, and a merge drops a crossing. ARCHITECTURE.md forbids a producer from self-certifying a mark. With `hits.rs:88` and `marks.rs:48`, these are six of the ten unsafe paths in [NOT-EXACT-INVENTORY.md](../NOT-EXACT-INVENTORY.md), section 2.2.
Migration: callers that accepted a marked curve-surface answer now receive the refusal. A producer with real evidence for a mark builds it with `Mark::new` and its own code, kind and bound.
