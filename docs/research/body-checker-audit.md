# Independent body checker audit

Outcome: audit done, checker not built. Closed by the owner's certification choice of 14 November.
When: 11 to 14 November 2026.

## What we were trying to achieve

Know what the promise "never wrong when it says proven" really rests on. A finished body passes a topology check. The question was how much of the geometry that check re-derives for itself, and what a checker that trusts nothing would cost.

## Why it was proposed

It is step R10 ([archived plan, section 12.6](../archive/PLAN-2026-11.md#126-robustness-and-many-core-steps)). The reason on record is one sentence: the promise is only as good as the check behind it. Only the first step, the audit, was adopted. Building the checker was "named for later", to be decided from the audit's numbers.

## What was tried

On the target parts P0 to P5:

- a count of which facts the topology check reads and who recorded them;
- a rebuild with every recorded fact removed, to see what a "facts-blind" checker would have to recompute and at what cost;
- four planted faults: an edge moved off its faces, a face turned inside out, two overlapping lumps and two coinciding spheres.

It was run once on 11 November and again at the exit of the tolerance milestone.

## Result

- The topology check re-derives none of the geometry. 816 of its 872 questions (94%) are answered by a fact that a builder recorded.
- All four planted faults pass the topology check.
- The facts-blind checker catches all four.
- Its cost ranges from 0.1 times the work of the build (P1) to 25 times (crossing pipes, 99.3M units). So checking can cost more than building.
- It does not finish on P5: 20 curve-on-surface questions refuse on the spline faces.

The audit recommended building the checker in three stages, to run in CI and on request, never on a click.

## Why it was closed

The audit said "build it". The owner's later choice went the other way. On 14 November, choosing between a general independent certificate checker and "report what the checks establish", he took the second. `VISION.md` now says Catallax does not build independent certificate systems.

So the recommendation was not taken, for a product reason and not because the audit was wrong. The planted faults remain a real finding: no independent check guards against that kind of fault today.

The record does not show a separate ruling that names R10. This roadmap reads the 14 November choice as closing it. If the first stage (edges lie on their faces, under 50,000 units per analytic body) is wanted as a CI test, it can be asked for at a checkpoint.

## Sources

- [R10-BODY-CHECKER-AUDIT.md](../R10-BODY-CHECKER-AUDIT.md): https://github.com/BohnBawerick/catallax/pull/442
- [milestones/M12.md](../milestones/M12.md), "R10 audit at the exit source"
- `VISION.md`, section "Efficiency, not machinery"
- The options report "Twelve product choices for Catallax", 14 November 2026, kept outside this repository
