# Threads inside one field request

Outcome: adopted. The field kernel uses every processor by default.
When: 11 to 12 November 2026.

## What we were trying to achieve

One field request (a display mesh, a print mesh, a grid for analysis) that uses all the processors and still gives the same bytes, the same work count and the same refusals at any thread count.

## Why it was proposed

It is step R4 of the robustness and many-core steps the owner adopted on 11 November ([archived plan, section 12.6](../archive/PLAN-2026-11.md#126-robustness-and-many-core-steps)). The reason on record: the field speed that analysis needs belongs to the field track, and the analysis grid and the print cells both come after this step. The same day the owner ruled that Catallax uses every core and the GPU by design, without gating that on speed studies.

An earlier probe had reported identical bits at 1, 2, 4 and 8 threads with a 5.7 to 6.0 times speedup.

## What was tried

Two designs for sharing one work allowance between threads:

- A shared counter that every thread draws from. Rejected before measuring: the schedule would then decide which input runs out of work first, so the answer would depend on timing.
- A shared worker pool with ordered admission. Threads compute in parallel. The coordinator splits the allowance into fixed shares per input, charges the work and appends results in the original input order. This is the one that was built.

One measured round: P1 to P4, display and print, at 1, 2, 4, 8 and 12 processors, on the unchanged 1M allowance.

## Result

50 case-runs gave identical output bytes, work counts, cell counts and refusals at every processor count. No tolerance or allowance was raised.

In that round the 12-processor total was slower than the single-processor one. The machine was shared with other jobs, so this is not a timing result.

## Why it was adopted

The promotion rule was correctness only, and it passed. Speed was not a gate, in line with the 11 November ruling. Timing on a quiet machine belongs to the speed phase.

## Sources

- [prototypes/field-threads.md](../prototypes/field-threads.md)
- https://github.com/BohnBawerick/catallax/pull/413
- [archived plan](../archive/PLAN-2026-11.md) section 12.6, row R4
