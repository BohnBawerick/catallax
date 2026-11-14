# Fillet prototypes: which kernel builds the fillet, and patch-first construction

Outcome: both rounds parked. No fillet method came out of them.
When: 22 to 27 October 2026.

## What we were trying to achieve

A general fillet the user can see within the 1,000 ms feature budget, with a face proved within 100 nm of its rolling-ball definition, names that survive a radius edit and a STEP file that opens in OpenCASCADE without repair.

The open question was which kernel should build the fillet the user sees: the B-rep, or the field with the B-rep following later for export.

## Why it was proposed

The plan of 22 October allowed a prototype only when the answer is unknown, the answer changes what gets built, and a throwaway version is far cheaper than a proper build ([archived plan, section 8.1](../archive/PLAN-2026-11.md#81-when-a-prototype-is-worth-it)). Fillet ownership met all three. It was prototype 1 in [section 8.2](../archive/PLAN-2026-11.md#82-prototypes-now), and the general fillet build was to wait for its verdict.

## What was tried

Round 1, fillet ownership. Four fixed cases: crossed cylinders (G1), an oblique cone cut (G2), a sphere and cylinder edge (G3) and the P3 boss as an exact reference.

- B-rep first: cover the exact centre path with stations, build an envelope of exact rolling-ball spheres with ordinary Booleans, then trim.
- Field first: show the rolling-ball contact strips through the field, and build the B-rep only when STEP is asked for.

Round 2, patch-first. Replace the sphere unions with contact rails, exact circular sections and a small fitted spline cover. Two layouts: 32 cubic patches (A) and 16 degree-seven patches (B).

## Result

Round 1:

- B-rep first used the whole 100,000,000-unit allowance on G1 and produced no body, no display and no STEP. Twelve envelope unions took 97.7% of the measured time.
- Field first meshed G1 in a median of 323 ms and 520 ms at the two radii, inside the budget. Its bound was 7,500 nm for the display triangles, which is not the 100 nm whole-face proof. G2, G3 and P3 were outside what it could build. It gave no blend names and its deferred STEP ran out of work.

Round 2:

- Layout B had sampled errors of 30 and 37 nm. Layout A had 1,384 and 1,629 nm.
- Both produced open patch sets with zero solids. There were no trimmed neighbour walls, no shared edges with the neighbours and no whole-face proof.

## Why it was parked

The decision rule written before the runs asked for a complete solid, a whole-face proof, edit identity and a valid STEP. Neither round met it. Both notes say the result does not show that these fillets are impossible or too slow, only that these constructions did not finish.

The 10 November decisions replaced the question instead of resuming it. Whichever kernel is faster shows a feature first, the B-rep owns export and drawings, and general fillets come from the certified families plus the spline-face blends of milestone 1.5 ([archived plan, section 11.5](../archive/PLAN-2026-11.md#115-items-superseded-on-10-and-11-november)). The plan's own words are that the prototype "decided nothing".

The record does not say why a third construction was not tried. It records only that the family-by-family fillet work took over ([fillet-coverage.md](fillet-coverage.md)).

## Sources

- [prototypes/fillet-ownership.md](../prototypes/fillet-ownership.md): protocol https://github.com/BohnBawerick/catallax/pull/28, verdict https://github.com/BohnBawerick/catallax/pull/40
- [prototypes/fillet-patch-first.md](../prototypes/fillet-patch-first.md): https://github.com/BohnBawerick/catallax/pull/59
- [milestones/M1.md](../milestones/M1.md), "Parked work and reasons"
- [archived plan](../archive/PLAN-2026-11.md) sections 8.1, 8.2 and 11.5
