# Tolerance as a speed lever

Outcome: dropped. The 1 nm rule stays for correctness, not for speed.
When: 10 November 2026.

## What we were trying to achieve

Find out where the work goes in the three P0 cases that run out of the 1,000,000-unit scoreboard allowance: block and sphere, crossing pipes and the seam case.

## Why it was proposed

The plan of 22 October said the expensive step is proving that something is exactly zero, and that any tolerance removes that step. If that were true, the 1 nm decision rule would be the way to make the slow cases fit. Nobody had measured it. The same change also produced the inventory of every inexact decision in the code, which the architecture required before the tolerance rules could be called enforced.

## What was tried

A profiling switch in the numbers crate that attributes every work unit to the decision that spent it, split by outcome. It reads the counter and changes no answer. The three cases were run at 1M and at 100M.

## Result

Exact-zero proofs are a small part of the work:

| Case | Share of the whole request |
|---|---:|
| Block and sphere | 12.4% |
| Crossing pipes | 0.4% |
| Seam | 2.1% |

With every exact-zero proof removed, and every enclosure that failed to separate zero, the three requests would still need 2.35M, 9.10M and 6.30M units. That is 2.4 to 9.1 times the 1M allowance.

The real costs are high-precision coefficient enclosures in the imprint probe path, root isolation of one intersection polynomial, and cancellation checks in the arctangent series during STEP export.

## Why it was dropped

The measurement showed the premise was false for these cases. The plan now says so ([archived plan, section 11.5](../archive/PLAN-2026-11.md#115-items-superseded-on-10-and-11-november)): the 1 nm rule stays for correctness, and the three cases move to refusal repair as "growing complexity" and "slowly converging evaluation".

Limit of the evidence: three P0 cases. The record has no such profile for spline parts.

## Sources

- [NOT-EXACT-INVENTORY.md](../NOT-EXACT-INVENTORY.md), section 1
- https://github.com/BohnBawerick/catallax/pull/359
- [archived plan](../archive/PLAN-2026-11.md) sections 3.1 and 11.5
- Related: [decide-once-merging.md](decide-once-merging.md)
