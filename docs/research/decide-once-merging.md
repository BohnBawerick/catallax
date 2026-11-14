# Decide-once merging

Outcome: dropped. Step 2 of the tolerance milestone (old M12) closed with no merging.
When: 11 November 2026.

## What we were trying to achieve

Fewer refusals from decisions that are too close to call. The first consistency rule of the plan is "decide once, everyone obeys": a too-close decision is recorded, and every later decision reads the record instead of asking again ([archived plan, section 3.3](../archive/PLAN-2026-11.md#33-four-consistency-rules)). Merging two things that are within 1 nm is one way to make such a decision.

## Why it was proposed

The tolerance milestone allowed merging only at the 22 places the inexact-decision inventory lists as candidates, and only where it changes an answer for the better. So the question was whether any saved answer would improve.

## What was tried

The 477 saved questions were run at their 100M allowance. Each result counts its undecided decisions. The 17 that ran out of work were profiled again to see where the work went.

## Result

None qualifies. No saved answer makes an undecided decision at any of the 22 places. Of the 17 that run out of work, three are dominated by exact-zero proofs, and those spend the work in ordering and orientation tests, which a tolerance cannot answer.

## Why it was dropped

Merging would have changed no answer, so it was not built. The three cases went to refusal repair. The decide-once record itself stays in the design for correctness.

## Sources

- [NOT-EXACT-INVENTORY.md](../NOT-EXACT-INVENTORY.md), section 4
- https://github.com/BohnBawerick/catallax/pull/404
- [archived plan](../archive/PLAN-2026-11.md) section 13.1, row M12
- Related: [tolerance-as-speed-lever.md](tolerance-as-speed-lever.md)
