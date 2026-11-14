# Census of 3,000 public STEP files

Outcome: recorded. It is the measure of the import gap.
When: 11 to 12 November 2026.

## What we were trying to achieve

Know what share of real parts Catallax can open today, and why the rest fail, sorted by cause.

## Why it was proposed

It is step R8 ([archived plan, section 12.6](../archive/PLAN-2026-11.md#126-robustness-and-many-core-steps)). Published numbers count faces, not parts: planes are about 61% of faces in one large sample, cylinders about 29%. None of them says what share of whole parts Catallax can carry today. The plan wanted that number to tell the owner how far the everyday-part milestones reach, since the target parts P1 to P10 were designed to exercise features and were not taken from anyone's work.

## What was tried

A first attempt on 11 November stopped before measuring, because the archive could not be fetched. The run that counted took 3,000 models from the public ABC collection with a fixed seed and put each through the import entry on main. It was run again after a fix to how placements are read.

## Result

| | Files | Share of 3,000 |
|---|---:|---:|
| Imported with every source item covered | 186 | 6.2% |
| Imported, some items not verified | 3 | 0.1% |
| Refused with a named cause | 2,782 | 92.7% |
| Timed out | 29 | 1.0% |

Split by what the file contains:

- 2,003 files have only planes, cylinders, cones, spheres and tori. 186 of them import in full.
- 868 files have spline faces. None of them imported in this run.
- 129 files could not be classified.

The most common first causes of refusal: a face the importer finds self-intersecting (822), an edge on a spline solid it cannot certify (532), a face on a spline solid it cannot certify (196), and exact incidence checks (203).

The placement fix removed 982 diagnostics and admitted no new file. Those files went on to fail a later check. The census stops at the first refusal, so it cannot say how many files a single fix would admit.

A follow-up on 12 November took the 728 files stopped by spline certification and removed that blanket gate:

- 412 then fail because a curve or point in the file does not lie exactly on its face.
- Every one of those 412 gaps was measured. All are below 0.006156 mm, under the 0.01 mm repair limit that import already has.
- One file now imports, as four marked bodies. The other 727 still refuse.

## What it changed

The note's own conclusion: spline solid certification, reported intersections, exact incidence and incomplete proofs are the main import boundaries.

The follow-up shows the pattern behind the largest group: the files are close, by a few micrometres, and an exact rule rejects them. That is why real STEP import is the first gap on this roadmap, and why it is built on numerical modelling under stated tolerances (milestones 1.1 and 1.2). The census is rerun at every checkpoint.

Not measured: fillet and chamfer counts (STEP carries no feature history), and whether each source file is itself valid.

## Sources

- [milestones/R8-public-step-census.md](../milestones/R8-public-step-census.md) and its table `R8-public-step-census.csv`: https://github.com/BohnBawerick/catallax/pull/415 and https://github.com/BohnBawerick/catallax/pull/430
- [milestones/R8-public-step-preflight.md](../milestones/R8-public-step-preflight.md)
- [milestones/R8-incidence-gaps.md](../milestones/R8-incidence-gaps.md) and [milestones/R8-spline-certification.md](../milestones/R8-spline-certification.md), the 12 November follow-ups
- The placement fix: https://github.com/BohnBawerick/catallax/pull/425
