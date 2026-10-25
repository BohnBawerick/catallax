## 2026-10-20: zero radicands retain working precision

Piece: 1
Change: `enclose_at` tries a bounded exact tower proof when a square-root argument's ball overlaps zero. A proved zero is enclosed as zero. Other known-defined arguments still use the square root of the ball's nonnegative part, without annihilator or separation-bound work.
Why: clipping an uncertain enclosure of an exact zero loses half the working bits at each nested square root. In the smoke probe, a nested zero root at 128 bits had width `3.1798e-10`; at 256 bits it still had width `5.6398e-20`. The exact proof gives width zero at both precisions and changes no decision statistics.
Migration: no signature, accuracy, refusal code or budget default changes. Positive radicands must never be clamped to zero solely because their balls overlap zero.

