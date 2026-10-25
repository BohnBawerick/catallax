## 2026-10-19: what a boolean stage record's statistics mean

Pieces: 5, cli
Change: in `StageRecord::stats`, every counter other than `max_bits` is what the stage added to the `Decider`'s statistics, and `max_bits` is the `Decider`'s running `max_bits` when the stage ended, an upper bound on the stage's own peak. Piece 5 reads statistics only through `stats()`, never `take_stats()`.
Why: piece 1's `Decider` keeps only a running maximum precision and offers no way to restore cleared statistics, so a stage's own `max_bits` cannot be observed without clearing the caller's totals.
Migration: piece 5 snapshots `dec.stats()` around each stage. Its test writer checks case 53 as now written. No signature changes.

