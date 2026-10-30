## 2026-10-30: P5 construction work inside the application allowance

Pieces: shapes, features, judge interface
Change: `RollingEnvelope::fit` now keeps each span coordinate as `even + t * odd` for the transition tangent `t`, so controls with rational node data no longer carry a deep expression over `t`. Full P5 construction drops from 24,540,684 to 2,002,748 work, and the radius 5/2 edit from 24,526,585 to 2,001,242. At 20M the ordinary P5 request therefore no longer returns `numbers.operation-budget`. It reaches the existing `numbers.geometric-gap-unproven` exchange refusal, as it already did at 100M.
Why: construction alone exceeded the 20M application allowance, before any exchange, export, display or measurement work.
Migration: a caller that treated budget exhaustion as P5's 20M outcome now sees the exchange refusal, with no STEP file. Fit nodes, spans and proved bounds are unchanged, and the carrier controls are the same exact values under smaller expressions. The 1M allowance still ends in `numbers.operation-budget`. No format, recipe or allowance changes.
