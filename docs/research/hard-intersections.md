# Hard-intersection prototype

Outcome: both candidates parked.
When: 22 October 2026.

## What we were trying to achieve

Booleans that finish on the hard analytic poses: equal-diameter cylinders crossing at right angles, cylinders tangent along a line, a cylinder and a sphere of equal radius, and poses along a cylinder seam. At the time these refused because the exact algebra grew into very large numbers.

## Why it was proposed

It was prototype 2 in [archived plan, section 8.2](../archive/PLAN-2026-11.md#82-prototypes-now). The plan named two ways to avoid the large polynomials and did not know which would work, so it asked for one measured round before the repair work started.

## What was tried

Two candidates against the unchanged main branch, on 56 fixed inputs (24 canonical cases and 32 requests from the hard-case survey), each with 100,000,000 work units and a 30 second guard.

- Implicit: keep the meeting curve as "where these two surfaces meet" and decide against that directly, without expanding the elimination polynomial.
- Certified approximation: keep the definition beside chord tubes with a proved error of at most 100 nm, and use the tubes only to prove that things do not meet.

## Result

| | Baseline | Implicit | Certified |
|---|---:|---:|---:|
| Answered, of 56 | 15 | 10 | 7 |
| Incomplete at the guard | 38 | 44 | 49 |
| Cases with a paired median inside the 500 ms Boolean budget | | 0 | 0 |

Neither candidate answered more cases than the code it was meant to replace. On the 477 earlier requests, implicit kept 148 earlier answers and left 29 incomplete. Certified kept 70 and left 107 incomplete. No completed answer changed.

## Why it was parked

The rule was to promote the candidate that answers the most cases inside the click budget with no changed answer. Neither did. The note keeps both branches parked and says the general and singular cases need a different construction or proof.

On 10 November the plan marked the prototype obsolete for spline faces, because candidates there now come from the ported SMLib intersection code (now milestone 1.3). The analytic hard cases went to refusal repair ([archived plan, section 11.5](../archive/PLAN-2026-11.md#115-items-superseded-on-10-and-11-november)).

## Sources

- [prototypes/hard-intersections.md](../prototypes/hard-intersections.md): protocol https://github.com/BohnBawerick/catallax/pull/24, verdict https://github.com/BohnBawerick/catallax/pull/27
- [milestones/M1.md](../milestones/M1.md), "Parked work and reasons"
- [archived plan](../archive/PLAN-2026-11.md) sections 8.2 and 11.5
