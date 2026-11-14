# Spline against spline: the 2026 method with a proof per cell

Outcome: adopted, with an exact proof under every cell.
When: 11 to 12 November 2026.

## What we were trying to achieve

Find every curve where two spline faces meet, small loops included, and know that none was missed. With no equation on either side the search runs in four unknowns. A missed loop gives a wrong solid.

## Why it was proposed

It is step R13 ([archived plan, section 12.6](../archive/PLAN-2026-11.md#126-robustness-and-many-core-steps)). Two recent papers claim a guaranteed search: Yang, Jia and Yan (2023) and Yang and Jia (2026). Both work in floating point. The open question was which of them Catallax can run in exact arithmetic inside the allowance, beside the older rule that no closed loop exists where the two normal cones do not overlap.

## What was tried

An experiment on 24 pairs: six kinds of pair (a wavy pair, a loop, a tangent contact, a small loop, a disjoint pair and an open control), each as a base, moved, turned and scaled copy. Three tests were run in exact arithmetic: normal-cone exclusion, the 2023 critical projection and the 2026 field.

## Result

The experiment's own verdict was cautious. The 2026 paper limits its guarantee to the spacing of its grid and says it can miss loops when several critical points share one cell. Exact arithmetic does not close that gap. So the experiment recommended the 2026 method only as a source of starting points, and named the 2023 route as the next thing to prototype.

## The decision

The owner decided on 11 November, after the three choices were laid out: use the 2026 method as the engine, with an exact proof layer under every cell. This replaced the experiment's "2023 route next".

What that means in the code: subdivision replaces the fixed grid. A cell is accepted only with an exact proof that it holds no unfound branch. Any other cell splits. A cell still unresolved at the depth or work limit refuses and names itself.

## What the build then showed

Every one of the 24 pairs gets its result inside the unchanged 20M allowance. The largest costs 16.9M units.

- 20 pairs are proved: branches found, or proved empty.
- The 4 tangent-contact copies refuse and name the exact contact point.
- A tiny loop about 0.000024 wide in parameter space is found in every pose. The paper's grid and the experiment's grid were both thousands of times too coarse to see it.

Still refusing at that point: a branch that runs along a face border, and a contact on a border.

## Sources

- [R13-SPLINE-TOPOLOGY.md](../R13-SPLINE-TOPOLOGY.md): https://github.com/BohnBawerick/catallax/pull/417
- [M10-SPLINE-SPLINE.md](../M10-SPLINE-SPLINE.md): https://github.com/BohnBawerick/catallax/pull/468
- Yang and Jia, ACM Transactions on Graphics 45(5), 2026: https://doi.org/10.1145/3807948
- Yang, Jia and Yan, ACM Transactions on Graphics 42(6), 2023: https://doi.org/10.1145/3618349
