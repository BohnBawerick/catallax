# Scan fitting: from points to a clean part

Outcome: first cut built. Learned segmentation not adopted.
When: 12 to 13 November 2026.

## What we were trying to achieve

Turn a scan, a cloud of measured points, into a clean workable part with the smallest deviation, as a B-rep and as a field, with no feature tree needed as input.

## Why it was proposed

On 10 November the owner ruled that both kernels must act on the same geometry, imported parts and scans included, and that a scan should become a clean part without a recipe. That replaced the earlier rule "never infer a solid from a scan". The owner also set the order of methods: classic segmentation and fitting first, a learned method such as DualBrep only if the classic route is not good enough or does not scale. Scan to part is now part of milestone 4.4, and its test is that 30,000 noisy points taken from P3 come back as P3's planes and cylinders.

## What was tried

Classic fitting. Estimate a normal at each point, grow flat regions, merge regions that one plane, cylinder or cone explains, then assemble: a block as stock, bosses added, bores cut, rounds on the stock edges. The output is an ordinary recipe that both kernels build. Deviation is reported both ways, scan to part and part to scan.

A research report then asked whether learned segmentation should be added.

## Result

The P3 round trip, 30,000 points with 0.01 mm noise:

| P3 surface | True | Fitted |
|---|---|---|
| Block faces | 0 and 40, 0 and 30, 0 and 10 | all within 0.0013 mm |
| Boss wall | radius 8 | radius 8.0005 |
| Top-edge rounds | radius 1 | radius 1.0059 to 1.0127 |
| Round where the boss meets the block | a torus, radius 2 | refused: 885 points left unexplained |

So it meets the milestone's words, planes and cylinders, and misses the one torus.

Answers and work counts are the same at 1, 2 and 4 workers. The deviation numbers are maxima over samples. They are not proved bounds.

## The verdict on learned segmentation

Do not add it now. The P3 miss is a missing torus fit and a missing assembly rule, not a segmentation failure: a perfect labelling of the points would still not let this assembler build a pocket or a chamfer. DualBrep is worth one bounded comparison on free-form cases later. No learned model was run, so the verdict is from reading, not from a measured gain or loss.

## What is still open

Tori, chamfers, pockets, steps and features off the stock axes. Free-form regions wait for spline Booleans and spline fillets (milestones 1.3 and 1.5). The record shows no run on a real scan, only on points sampled from our own parts.

## Sources

- [M23-SCAN-FIT.md](../M23-SCAN-FIT.md): https://github.com/BohnBawerick/catallax/pull/502 and https://github.com/BohnBawerick/catallax/pull/503
- [archived plan](../archive/PLAN-2026-11.md) section 11.5, the scan row
- The research report "M23 learned-segmentation verdict", 13 November 2026, kept outside this repository
