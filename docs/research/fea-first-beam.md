# First FEA beam

Outcome: parked by the owner. The cantilever runs. Its numbers are unproved.
When: 11 to 13 November 2026.

## What we were trying to achieve

Stress and deflection computed on the same part the B-rep exports, using the field kernel's grid. FEA is the stated goal of the field track.

## Why it was proposed

On 11 November the owner asked for the simplest possible case first: a rectangular cantilever beam, one face fixed, a load on another face, then the von Mises stress along the beam and its deflection. In his words, getting that done with the field would be a good start.

The plan treats FEA as a proof of the pipeline, with a design note written before the groundwork finished. The design had one known risk to cover: grid methods are weakest at the boundary, where a cell cut small by the surface degrades the stiffness matrix.

## What was tried

1. A first cut: three-dimensional linear elasticity on a plain grid, on a 100 by 10 by 10 mm beam.
2. A design note for the first field cantilever and the later stages.
3. The agreed cantilever from that design, with stress recovery, a stress picture on a slice and the centreline deflection drawn beside beam theory.

## Result

The first cut's deflection is within 0.89% of beam theory at 8 by 8 by 80 cells. A coarser 4 by 4 by 40 grid through the analysis adapter is within 3.5%.

Every displacement, stress and reaction is an unproved numerical estimate. Nothing here feeds a proven answer or an export. The plate, the sliver case and an imported part are separate acceptance lines of the milestone, and this work does not claim them.

## Why it was parked

The owner parked all FEA work on 13 November: finish only what was close to done, start nothing new. The agreed cantilever was close and merged that day. The stress study in the app and two FEA feature ideas were parked with it.

The record holds the ruling but not a reason in the owner's words. It came the same day as the push to close the release, and this roadmap files it under "finish modelling first". FEA restarts on the owner's word at a checkpoint.

## Sources

- [FEA-FIRST-CUT.md](../FEA-FIRST-CUT.md): https://github.com/BohnBawerick/catallax/pull/390
- [FEA-DESIGN.md](../FEA-DESIGN.md): https://github.com/BohnBawerick/catallax/pull/385
- [M20-BEAM.md](../M20-BEAM.md): https://github.com/BohnBawerick/catallax/pull/568
- [archived plan](../archive/PLAN-2026-11.md) section 12.6, row R14, on cut cells
