# Cheaper proof on P1

Outcome: adopted. P1 fits the 20M allowance with no change to any rule.
When: 13 November 2026.

## What we were trying to achieve

P1, the bracket, building inside the app's unchanged 20,000,000-unit allowance. It needed 41.6M. The wider aim was to learn which kinds of proof saving are worth building.

## Why it was proposed

The plan named "cheaper proof" as research: proof from recorded construction facts, reuse of proofs, proof scoped to what was asked for, and coarse proofs refined where needed ([archived plan, section 8.2.1](../archive/PLAN-2026-11.md#821-cheaper-proof-prioritised-research)). It was the only planned help for the proof half of the click targets.

On 13 November the owner also ruled that kernel work must be general: a change must apply to every part, and tuning one sample to fit under 20M is worthless. The record does not say the P1 change was started because of that ruling. It is the standard the change was read against in the proof-cost research the same day ([proof-cost.md](proof-cost.md)).

## What was tried

Three changes, measured one after the other on the eager build of P1:

| Change | Work after it | Saved |
|---|---:|---:|
| Starting point | 41,647,298 | |
| An exact proof for parallel cylinders, from how they were constructed | 33,754,179 | 7,893,119 |
| Reuse of one separating enclosure | 29,866,897 | 3,887,282 |
| Reuse of enclosures inside one request | 15,264,436 | 14,602,461 |

None of them swaps a proof for a sample. The construction proof keeps its stored witness and proves clearance along the whole path. The reuse keeps each decision's own arithmetic limit.

## Result

41,647,298 to 15,264,436 units, a saving of 63%. P1 has 4.7M units of headroom under 20M. Most of what remains is in intersection (54%) and imprint (28%).

It does not meet the 2,000 ms proof target and does not complete P1's acceptance in the app. Those are separate gates.

## Why it was adopted

It met the allowance by doing less repeated work, with every guarantee kept, and the methods are general. It merged on 13 November.

A research report written beside it ranks what to build next, and this roadmap keeps that order in the speed phase (milestone 5.1): proofs from recorded construction facts, prepared polynomial bounds, reuse inside a request, then reuse across edits for regions that did not change, then coarse proofs refined only where they fail.

## Sources

- [M2-P1-WORK-PROFILE.md](../M2-P1-WORK-PROFILE.md): https://github.com/BohnBawerick/catallax/pull/546
- [archived plan](../archive/PLAN-2026-11.md) section 8.2.1
- The research report "M26 cheaper proof: measured costs and build recommendation", 13 November 2026. It is kept outside this repository.
