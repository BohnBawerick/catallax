# Adaptive printing of field parts

Outcome: parked in parts. It works on a convex block. The proof it needs for cavities does not exist yet.
When: 11 to 13 November 2026.

## What we were trying to achieve

A printable mesh of a lattice part at the default 0.01 mm accuracy. The test part is P6, a small bracket with a 1 mm shell and 2 mm lattice cells inside.

## Why it was proposed

A uniform grid cannot do it. At 0.01 mm, one layer of cells on a single 40 mm² face of P6 already needs about 10.5 million cells, against a limit of 2,000,000. The standing print rule says field parts use adaptive cells: small cells only where the surface needs them. Milestone 4.3, lattices and printing, is the milestone for that.

## What was tried

An opt-in adaptive print request: cells that refine where needed and stay conforming, two independent distance checks on the result, and up to four rounds (one first try and three refinements). Exact point queries for sheet shapes such as gyroids were added with it.

## Result

It is admitted for one kind of shape with material in it: a single convex block. Every other shape refuses at the last check, the one that asks whether the mesh has the same pieces, cavities and handles as the definition.

That refusal is correct. The verdict note gives an exact counterexample: a block with a spherical void inside. A mesh can be closed and within the distance bound on both sides and still have lost the void. Distance alone does not prove the material has the right shape.

The other parts of the lattice and printing milestone at the verdict:

- Sheet gyroid printing: parked until the material proof exists.
- P6 at 0.01 mm: does not fit the current limits on the uniform route.
- Cost estimates for a coarser choice: parked. There are no calibration samples yet.

The older P6 print at 0.5 mm still works.

## Why it was parked in parts

The release rule let research lines finish with a written verdict. The guard that refuses could not be removed on the evidence, because the counterexample shows it would let wrong meshes through. What was built stays in the kernel: the adaptive request, the sheet queries and failure reports that say how far a run got.

What is missing is a proof that a candidate mesh has the same material layout as the definition. That is the first job of milestone 4.3 when the field phase starts.

## Sources

- [parked/m19-printing.md](../parked/m19-printing.md): https://github.com/BohnBawerick/catallax/pull/575
- The adaptive increment: https://github.com/BohnBawerick/catallax/pull/405
- [archived plan](../archive/PLAN-2026-11.md) section 3.9, print accuracy
