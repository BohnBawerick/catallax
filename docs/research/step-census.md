# Census of 3,000 public STEP files

Outcome: recorded. It is the measure of the import gap.
When: 11 to 15 November 2026.

## What we were trying to achieve

Know what share of real parts Catallax can open today, and why the rest fail, sorted by cause.

## Why it was proposed

It is step R8 ([archived plan, section 12.6](../archive/PLAN-2026-11.md#126-robustness-and-many-core-steps)). Published numbers count faces, not parts: planes are about 61% of faces in one large sample, cylinders about 29%. None of them says what share of whole parts Catallax can carry today. The plan wanted that number to tell the owner how far the everyday-part milestones reach, since the target parts P1 to P10 were designed to exercise features and were not taken from anyone's work.

## What was tried

A first attempt on 11 November stopped before measuring, because the archive could not be fetched. The run that counted took 3,000 models from the public ABC collection with a fixed seed and put each through the import entry on main. It was run again after a fix to how placements are read.

## Result

The complete 15 November rerun imports **276 of the same 3,000 files with full
source coverage: 9.20%, beside the previous 6.20% (about 6%)**. All 3,000 inputs
were observed, with the original seed, source-byte hashes and denominator.

| Outcome | 11 November files | Share | 15 November files | Share |
|---|---:|---:|---:|---:|
| Imported with every source item covered | 186 | 6.20% | 276 | 9.20% |
| Imported, some source items not verified | 3 | 0.10% | 3 | 0.10% |
| Refused with a named cause | 2,782 | 92.73% | 2,719 | 90.63% |
| External observation timeout | 29 | 0.97% | 2 | 0.07% |
| Total observed | 3,000 | 100% | 3,000 | 100% |

All 186 previous full imports and all three partial imports remain in those
categories. Ninety previously refused files now import fully; none of the old
timeout cases contributes an admission gain. Twenty-seven old timeouts now
return a refusal and two remain timed out. The returned categories are 190
imports through the existing certified path and 89 numerical approximate imports;
three of the certified imports retain unverified source coverage.

The [per-model table](step-census/results.csv) retains each original first cause
beside the current outcome and cause. The [run manifest](step-census/results.json)
pins the source revision, observer, collector, binary, toolchain, caps and table
hash. The ninety recovered full imports were previously stopped by
`self_intersecting` (83), `not_solid` (6), and spline-edge certification (1).
These are actual recovered model IDs, listed in the manifest.

### First causes before and after

| Returned cause | 11 November | 15 November |
|---|---:|---:|
| Planar source repair not verified | 0 | 975 |
| Source incidence gap | 0 | 505 |
| Unsupported source type or variant | 1,314 | 876 |
| Resource limit | 334 | 126 |
| `not_solid` | 203 | 80 |
| Self-intersection not verified | 85 | 79 |
| `self_intersecting` | 822 | 53 |
| Malformed input | 24 | 24 |
| Numerical family unsupported | 0 | 1 |
| Total returned refusals | 2,782 | 2,719 |

The manifest preserves every exact diagnostic bucket, including unsupported
types and variants. Changed first causes expose later checks: the reduction from
822 to 53 reported self-intersections corresponds to 83 actual admissions from
that original group, rather than 769. Missing repair proofs, source limits and
unsupported families remain explicit refusals.

All 89 numerical imports use the fixed 0.01 mm tolerance and remain marked
approximate. Their largest recorded local source-to-result bound is
`1/8000000000` mm. That source limit is separate from the numerical body's checks;
a complete-body approximation bound remains unestablished. This measures import
admission. The separate [oblique STEP qualification](../phase1/import-oblique-step/README.md)
checks actual pose and export workflows on fixed real sources.

### Earlier observations, 11 to 12 November

The 11 November run, split by what the file contains:

- 2,003 files have only planes, cylinders, cones, spheres and tori. 186 of them import in full.
- 868 files have spline faces. None of them imported in this run.
- 129 files could not be classified.

Its most common first causes of refusal were a face the importer found self-intersecting (822), an edge on a spline solid it could not certify (532), a face on a spline solid it could not certify (196), and exact incidence checks (203).

The placement fix removed 982 diagnostics and admitted no new file. Those files went on to fail a later check. The census stops at the first refusal, so it cannot say how many files a single fix would admit.

A follow-up on 12 November took the 728 files stopped by spline certification and removed that blanket gate:

- 412 then fail because a curve or point in the file does not lie exactly on its face.
- Every one of those 412 gaps was measured. All are below 0.006156 mm, under the 0.01 mm repair limit that import already has.
- One file now imports, as four marked bodies. The other 727 still refuse.

## What it changed

Spline solid certification, reported intersections, source incidence and incomplete proofs are the main import boundaries.

The follow-up shows the pattern behind the largest group: the files are close, by a few micrometres, and an exact rule rejects them. That is why real STEP import is the first gap on this roadmap, and why it is built on numerical modelling under stated tolerances (milestones 1.1 and 1.2). The census is rerun at every checkpoint.

Not measured: fillet and chamfer counts (STEP carries no feature history), and whether each source file is itself valid.

The SolidDraft half of milestone 1.2—replacing an imported file and re-pointing
every reference—remains open on the app side. This census is the coverage result.

## Input identity and reproduction

The seed remains `20261111`, sampling the original 3,000 ABC model IDs from
`abc_0000_step_v00.7z`. Its SHA-256 is
`52e6dd1b6fa38e3cd99af59b662370829129540030975919c3e256dce6ad1dbe`.
The [individual source hashes](step-census/source-hashes.json) have SHA-256
`3323e104c377023846073a05d7e164ca1082dc3dcf32cc046eee897401454cc7`.
All 3,000 hashes were independently checked after transfer to the work PC, and
each was checked before and after its observation. The
[ten fixed P06 capability sources](step-census/retained-controls.json) remain
retained privately with their original bytes and hashes; all ten import fully.
Raw geometry and archive member names stay private under the
[existing licence scope](../milestones/R8-public-step-census.md#source-and-licence-scope).

The observation source is
[`b83edf42b0dad259b738d1bf07a690f7c638464e`](https://github.com/BohnBawerick/catallax/tree/b83edf42b0dad259b738d1bf07a690f7c638464e),
whose production crate sources match main
`c0056400973a65cd4f9c4193a6436b4f157cdbf6`. The snapshot adds the observer's
opt-in numerical mode and collection tooling. This publication changes only
documentation and evidence. The comparison includes kernel work landed since
the 11 November baseline; the source revision identifies the complete observed
implementation.

From that pinned observation snapshot, build `catallax-import`'s `r8_census`
example in release mode and run:

```sh
python3 -B tools/research/import_census.py --binary TARGET/release/examples/r8_census --sources PRIVATE_SOURCES --selection ORIGINAL_PRIVATE_SELECTION.json --source-hashes docs/phase1/import-census/source-hashes.json --out NEW_PRIVATE_CAPTURE --workers 4
```

The release run uses Rust 1.100.0 on the Linux work PC, four independent observers
with one kernel worker each, the unchanged 20M per-import allowance, and the
unchanged external 60-second observation cap. `full` requires returned bodies and
complete source coverage; `partial` retains unverified source coverage. Timeouts
remain separate even if a checkpoint was emitted before termination. A complete
capture requires all 3,000 terminal observations; process errors or missing
results prevent publication as complete. Raw checkpoints remain private.

Validation includes the full Linux format gate, four observation controls
(partial coverage, external timeout, process failure and numerical category),
and the CI-selected eighteen-crate all-target Clippy gate on the observer source.
The published table was independently checked for its hash, 3,000 unique original
IDs, every pinned source hash, aggregate counts and retained controls.

## Sources

- [15 November per-model table](step-census/results.csv), [run manifest](step-census/results.json), [source hashes](step-census/source-hashes.json), and [retained controls](step-census/retained-controls.json)
- [milestones/R8-public-step-census.md](../milestones/R8-public-step-census.md) and its table `R8-public-step-census.csv`: https://github.com/BohnBawerick/catallax/pull/415 and https://github.com/BohnBawerick/catallax/pull/430
- [milestones/R8-public-step-preflight.md](../milestones/R8-public-step-preflight.md)
- [milestones/R8-incidence-gaps.md](../milestones/R8-incidence-gaps.md) and [milestones/R8-spline-certification.md](../milestones/R8-spline-certification.md), the 12 November follow-ups
- The placement fix: https://github.com/BohnBawerick/catallax/pull/425
