Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

## 2026-11-21: proven stationary pairs where no closed form holds the minimum

Pieces: measure (`catallax_measure`)
Change: `catallax_measure::solid_gap` settles the strata pairs that have no closed form (a line in general position against a circle, circles with non-parallel normals, and lines, circles and cones against cones) by their stationary pairs, wherever the certified interval search leaves their face pair unsettled. A Krawczyk test proves each stationary pair, and an exact pair of carrier points beside it is classified on both trims. Clear answers that stated a gap wider than the requested accuracy, with the `measure.accuracy-not-reached` mark, now state it within the accuracy and without the mark. No signature changes.
Why: roadmap 3.3, step 3. The interval search's lower bound converges only as fast as its cells shrink, so pairs with a cone, or tori at non-parallel axes, kept the strata boxes' separation as the lower end of their gap.
Migration: none. Every answer that was within the accuracy before keeps its enclosure and work. Answers that carried the mark may now give a narrower enclosure without it, at more work. A pass that runs out of allowance returns its work and leaves the step 1 answer.
