## 2026-10-21: certify stalled algebraic-zero decimals

Piece: 7
Change: after decimal refinement stalls on an enclosure containing zero, an exact piece 1 coincidence proof may establish the decimal `0.` with zero error. A containing interval alone never permits this result.
Why: a nonnegative enclosure of an algebraic zero can have a positive midpoint at every precision. Relative decimal refinement then stalls even though the defining polynomial proves the value exactly zero.
Migration: no signature, rounding rule or default budget changes. Previously successful certificates take the same path. Regressions use a cubic algebraic identity and a tiny positive perturbation to distinguish exact zero from a small unresolved value.


