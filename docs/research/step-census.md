# Census of 3,000 public STEP files

Outcome: recorded. It is the measure of the import gap.
When: 11 to 16 November 2026.

## What we were trying to achieve

Know what share of real parts Catallax can open today, and why the rest fail, sorted by cause.

## Why it was proposed

It is step R8 ([archived plan, section 12.6](../archive/PLAN-2026-11.md#126-robustness-and-many-core-steps)). Published numbers count faces, not parts: planes are about 61% of faces in one large sample, cylinders about 29%. None of them says what share of whole parts Catallax can carry today. The plan wanted that number to tell the owner how far the everyday-part milestones reach, since the target parts P1 to P10 were designed to exercise features and were not taken from anyone's work.

## What was tried

A first attempt on 11 November stopped before measuring, because the archive could not be fetched. The run that counted took 3,000 models from the public ABC collection with a fixed seed and put each through the import entry on main. It was run again after a fix to how placements are read.

## Result

The complete frozen-population export comparison records **369 of the original
3,000 files (12.30%) with STEP exports valid after OCCT builder construction,
which can heal pcurves internally**, beside **2,701 of 3,000 (90.0%) valid after
normal OCCT transfer and default ShapeFix healing**. Neither number establishes
unhealed validity, source fidelity or a common-tolerance comparison.

The conditioning step is **400 admitted files, of which 369 are valid after
builder construction (92.25%), 14 invalid, 15 export-refused and two timed out**.
All 400 imports reproduce their expected body counts (808 bodies).
A file counts only when every expected body exports and passes whole-solid
validity after builder processing.
The [complete comparison](step-census-circular-pcurves/README.md),
[3,000-ID table](step-census-circular-pcurves/results.csv),
[candidate manifest](step-census-circular-pcurves/results.json),
[matched-control manifest](step-census-circular-pcurves/matched-control.json) and
[paired outcomes](step-census-circular-pcurves/comparison.json) retain all observations.

The isolated constant-height circular-pcurve repair recovers **eight actual
files with no losses** from the matched integrated control's 361 valid files.
It preserves existing flat maps and tilted carriers, and retains the finite U
chart. All 699 written body-file 3D graphs, oriented loops and pcurve surface
associations remain identical; 33 byte streams change. Candidate and control
import work, source hashes, body counts and categories agree on all 400 inputs.
The only changed production path is the Plane/Cylinder numerical STEP writer.
A broader 363-file experiment was withheld after four gains and two losses.
The prior [seam-order comparison](step-census-cylinder-seams/README.md) retains
its distinct 354 to 361 measurement and seven verified recoveries.

All 15 returned refusals are `step.numerical-gap`. File `00002204`, previously
refused at the aggregate export work limit, now times out during export and
earns no recovery credit. `00004777` remains an export timeout. Both files'
export work and validity remain unknown.
All ten retained P06 controls remain valid after builder construction, including
its possible internal healing. Actual part `00009338` retains one solid,
matching the source-reference volume and all 45 material probes.
The independent audit checks all 400 durable source hashes and 699 declared
written body files; 677 reach builder checks, with 645 valid and 32 invalid bodies.
Valid bodies inside invalid multi-solid files earn no file credit.

The measured repair producer is
**`2d53b8b1bc307b075c6ff8323aed0595a9e93468`**, with executable SHA-256
`8b7e4f17a7ec3c0cedc3ba940d14c6fb9708a2b4de6e105d06399f45cc57212e`.
The matched control producer/kernel is
**`1e8b24b651c468e8f81cc1748fcc6a8b3fe3db3b`**. Its CSV and all 699 STEP files
reproduce the prior 361-file capture after the intervening main integration.
The prior [322-file export capture](step-census-export-validity/README.md) remains
unchanged with its distinct producer and protocol identity.
The original seed, source hashes and 3,000-file denominator are unchanged;
the 2,600 inputs outside the frozen admissions were not reevaluated for new
admissions. Later publication/rebase heads are distinct from measured producers.

The earlier **admission-only** census has 400 full source-coverage imports
(13.33%), beside the previous 378 (12.60%), 276 (9.20%) and original 186
(6.20%, about 6%). Its complete 3,000 terminal observations remain unchanged at
producer `51e9c0ffa13efcbcfb5c004b152e4c6014c6b5c1`. Historical admission counts
below do not establish post-builder export validity for those earlier producers.

| Admission outcome | 11 November | Prior checkpoint | Analytic repair | Periodic repair |
|---|---:|---:|---:|---:|
| Imported with every source item covered | 186 (6.20%) | 276 (9.20%) | 378 (12.60%) | 400 (13.33%) |
| Imported, some source items not verified | 3 (0.10%) | 3 (0.10%) | 3 (0.10%) | 3 (0.10%) |
| Refused with a named cause | 2,782 (92.73%) | 2,719 (90.63%) | 2,617 (87.23%) | 2,595 (86.50%) |
| External observation timeout | 29 (0.97%) | 2 (0.07%) | 2 (0.07%) | 2 (0.07%) |
| Total observed | 3,000 | 3,000 | 3,000 | 3,000 |

**All 378 previous full imports remain full; 22 actual additional files import
fully.** Every recovered file has recorded explicit cylinder seam repairs: 44
seams across these 22 files. Their original first refusals were 16
`self_intersecting`, three `not_solid` and three resource limits. Their actual
IDs, aggregate original causes and CSV digests are retained in the
[recovery summary](step-census-periodic/analytic-transition.json). Individual
source hashes are in the [original source manifest](step-census/source-hashes.json);
per-file outcomes and causes are in the [prior analytic table](step-census-analytic/results.csv)
and [periodic admission table](step-census-periodic/results.csv). The
[admission manifest](step-census-periodic/results.json) records all original/prior
outcome transition totals.
This is 124 additional full imports over 9.20%. The admission categories are 190
certified imports (including the three partial imports) and 213 numerical
approximate imports. Relative to 11 November, 213 refused files and one timed-out
file now import fully. All ten selected [P06 real controls](step-census-periodic/retained-real-controls.json)
remain full with unchanged source bytes.

The [admission per-model table](step-census-periodic/results.csv) retains each
original first refusal and the prior 9.20% outcome/cause beside the new outcome,
including each explicit topology repair and unresolved source-pair coverage.
The [admission manifest](step-census-periodic/results.json) binds the clean producer,
actual Cargo artifact, observer, collector, toolchain, caps and table hash.
The [12.60% table](step-census-analytic/results.csv) and
[manifest](step-census-analytic/results.json), and earlier
[9.20% table](step-census/results.csv) and [manifest](step-census/results.json),
remain unchanged.

The [initial analytic capture](step-census-analytic/initial-results.json)
[reported 382 full imports](step-census-analytic/initial-results.csv) (12.73%). After
restoring independent source-pair checks, four of those files returned
`self_intersecting`: `00001233`, `00005142`, `00007233`, `00007433`.
[Their identities and causes](step-census-analytic/initial-analytic-withdrawals.json)
remain retained. The historical 12.73% is not the current result. That initial
capture predates the actual Cargo-artifact binding required for current runs.
Intermediate incomplete runs and an unpublished periodic prototype census
contribute no current census claim.

### First causes before and after

| Returned cause | 11 November | Prior 9.20% | Analytic 12.60% | Periodic 13.33% |
|---|---:|---:|---:|---:|
| Planar source repair not verified | 0 | 975 | 129 | 129 |
| Source incidence gap | 0 | 505 | 502 | 591 |
| Unsupported source type or variant | 1,314 | 876 | 846 | 739 |
| Resource limit | 334 | 126 | 171 | 177 |
| `not_solid` | 203 | 80 | 80 | 80 |
| Self-intersection not verified | 85 | 79 | 78 | 80 |
| `self_intersecting` | 822 | 53 | 101 | 101 |
| Malformed input | 24 | 24 | 24 | 34 |
| Numerical family unsupported | 0 | 1 | 645 | 356 |
| Numerical finite trim rejected | 0 | 0 | 41 | 308 |
| Total returned refusals | 2,782 | 2,719 | 2,617 | 2,595 |

These are first stops, not predicted admission gains. Every exact unsupported
variant and transition remains in the manifests. This producer includes the
merged rational spline extrusion reader, so changes in unsupported/next-cause
buckets alone do not establish an import gain. The 22 recovered files are
identified independently by their actual admissions and seam-repair records.
Initial source bytes are charged before kernel verification hashing, including
mismatches. This can expose a resource limit before a later cause.

The finite analytic fallback retains shared 3D vertices/carriers and rebuilds
finite parameters. Its first periodic family joins two complete coaxial,
equal-radius circle rims with opposite winding and aligned existing anchors,
using one derived seam with two UV branches. Repeated cylindrical edges must
satisfy the full periodic-cell condition; a just-under-one-turn shared seam
cannot take the ordinary subperiod path. Complete circular holes require
continuous projected-disk/outer-segment clearance above the requested tolerance,
not only a sampled rim polygon. Unsupported multi-boundary and other periodic
families still refuse.

Cheap source checks run across every solid before optional curved proofs. The
latter share one collective 1M-work prepayment inside the unchanged 20M request.
Completed checks and every unresolved pair remain explicit. Neither unresolved
coverage nor zero sample residuals grants a global proof. Numerical STEP remains
approximate. Two follow-ups remain: less than 1M work remaining can cause refusal
at the collective prepayment, and legacy planar intersection work is not fully
metered. Counters cover existing reservations and added traversal/setup
checkpoints, rather than every legacy planar operation.

All 213 numerical imports use the fixed 0.01 mm tolerance. Their samples do not
establish global between-sample separation or a complete source/body bound. The
new circular-hole clearance is a separate bounded local check. Existing local
source-to-derivative bounds and pointwise residuals retain separate columns; an
absent source bound is not zero. These statements describe the admission capture.
The current export qualification adds structural post-builder STEP checks for its frozen
400 full inputs; it does not establish an arbitrary edit, material agreement or
a complete source/body bound for every file. The
[periodic qualification](../phase1/import-periodic-numerical/README.md) retains
the matching 51e9 named-part import/pose/STEP evidence and the later reviewed
two-seam repair, with their separate producer identities. Earlier
[analytic](../phase1/import-analytic-numerical/README.md) and
[oblique](../phase1/import-oblique-step/README.md) captures remain historical.

### Standard OCCT comparison on all 3,000 files

The complete reference run uses `STEPControl_Reader` normal root transfer,
`ShapeFix_Shape::Perform` with default settings, then `BRepCheck_Analyzer`, via
OCP 8.0.1.1. Source bytes are checked before and after every observation. The
[frozen protocol](../../tools/research/occt_census.py),
[reference manifest/default profiles](step-census-periodic/occt-reference.json),
[per-file comparison](step-census-periodic/comparison.csv), and
[comparison manifest](step-census-periodic/comparison.json) retain its scope.

| Exclusive reference outcome | Files | Share |
|---|---:|---:|
| Valid whole healed result, containing solids and no free lower-dimensional topology | 2,701 | 90.03% |
| Open-shell result with no solid | 44 | 1.47% |
| Other non-null opened result | 243 | 8.10% |
| Failed | 1 | 0.03% |
| External timeout | 11 | 0.37% |
| Total observed | 3,000 | 100% |

Opening a non-null transferred shape was observed for **2,993 files**, including
four that later timed out and one that failed during `shape_fix`. Thus 2,988
completed opened results plus four opened timeouts and one opened failure make
2,993 observed openings. An opening checkpoint does not establish completed
healing, validity, complete source coverage or a solid. The other seven timeouts
have no observed opening. Root transfer counts and whole-shape/solid checks are
retained separately. Empty cells mean unobserved, not a successful check.

All 400 Catallax full-import source files produce reference valid-solid results.
Another 2,299 files with a Catallax refusal and both Catallax timeouts produce
such results. The three Catallax partial imports produce other opened reference
results. All 44 reference open-shell cases remain Catallax refusals. These are
actual per-file outcomes, not a predicted conversion count.

The reference uses default reader/healing tolerances, not Catallax's 0.01 mm
pointwise cap, and carries no Catallax 20M accounting. It uses two file workers;
Catallax uses four file workers with one kernel worker each. Both retain a
60-second external per-file cap. This is a coverage comparison, not a speed or
common-tolerance benchmark. Default transfer already processes source geometry;
post-transfer validity is not raw source validity. Explicit ShapeFix changes one
observed invalid result to valid and five observed valid results to invalid;
unknown/timeout checks remain unknown.

For the owner's original 9.20%-checkpoint repair cohorts:

| Prior first refusal | Files | OCCT opening observed | Reference valid-solid result |
|---|---:|---:|---:|
| Planar repair not verified | 975 | 974 | 964 |
| Source incidence gap | 505 | 505 | 474 |
| `self_intersecting` | 53 | 53 | 53 |
| Self-intersection not verified | 79 | 79 | 78 |
| `not_solid` | 80 | 80 | 79 |

This provides actual files for subsequent healing work. Shapes that do not heal
into a closed solid require an explicit surface-body result; they are not counted
as Catallax full solid imports here.

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

Seed `20261111` selects the original 3,000 ABC model IDs from
`abc_0000_step_v00.7z`, SHA-256
`52e6dd1b6fa38e3cd99af59b662370829129540030975919c3e256dce6ad1dbe`.
The [individual source-hash manifest](step-census/source-hashes.json) has SHA-256
`3323e104c377023846073a05d7e164ca1082dc3dcf32cc046eee897401454cc7`;
the private original selection has SHA-256
`1b11410b220d5e8a3156794f6c8f9cc09091a43dcac9af801215fb3a38d571a5`.
All captures read the durable private source corpus. Every source hash is
checked before and after its observation. The
[ten fixed P06 capability sources](step-census/retained-controls.json) remain
retained privately with original bytes and hashes; all ten still import fully and pass the current post-builder export check.
Only public model IDs, hashes, diagnostic scopes and aggregate metadata are
published, under the [existing licence scope](../milestones/R8-public-step-census.md#source-and-licence-scope).

The measured full-admission producer is
[`51e9c0ffa13efcbcfb5c004b152e4c6014c6b5c1`](https://github.com/BohnBawerick/catallax/tree/51e9c0ffa13efcbcfb5c004b152e4c6014c6b5c1).
Its actual Cargo observer artifact SHA-256 is
`d6d645568442d5943e3325a209a98ef2b386077aca623de21a9ab63de6d1e6e2`.
Its admission per-file CSV has SHA-256
`b08224a6f919405de46d35896d8d682febbd10c36495d34759f4958c80be3bd7`.
The current export comparison binds candidate producer
`3651cd367cbb2f5d28118b06c6431193eb1b0d62`, actual executable
`17b928c6ee2bb18e522f3676a7db52fdfd64f2bb6dc6b323f137c87a3b2917bc`, and CSV
`0f9c1727baa480e4c4c9e446e9495fcfa891dfff0249a7acfbddbcc3e6928431`.
The matched control producer `a207ea676f52ee077cbb186674407dea3e64e1d0` has executable
`1c304cdcc090508bd8d9a533c37ffd6009263a0d862c7529e64e4b1019b4446f` and CSV
`8e765971a02100100e4c2a1bb7b81046bd9b2d01d763f7a1ab8b13baf318be3f`.
The separate import and aggregate export allowances remain 20M each under the
same 60-second whole-file observation.
The [comparison protocol and qualification](step-census-arc-parameters/README.md)
retain actual build identity, commands and work per stage.

The source/body fallback is numerical and opt-in at the kernel API. SolidDraft's
numerical-group and replacement-reference UI work remains with the app slice.

To reproduce the historical full-admission capture from its clean pinned producer, run:

```sh
python3 -B tools/research/build_import_observer.py --record PRIVATE_BUILD.json
python3 -B tools/research/import_census.py --binary ACTUAL_CARGO_ARTIFACT --build-record PRIVATE_BUILD.json --sources DURABLE_PRIVATE_SOURCES --selection ORIGINAL_PRIVATE_SELECTION.json --source-hashes docs/research/step-census/source-hashes.json --out NEW_PRIVATE_CAPTURE --workers 4
```

The builder records the actual compiler artifact and clean revision/tree, binary,
observer, lockfile and toolchain. The collector requires that identity. Its
`public-summary.json` is a whitelisted extract: local paths and full Cargo paths
remain private. `build_record_sha256` identifies the original private build
record; the nested public `build_record` contains its sanitized fields.

The Linux native release run uses Rust 1.99.1 (`58cf8f9c3`, 13 October 2026),
four independent file observers, one kernel worker each, the unchanged 20M
per-import allowance and the unchanged external 60-second cap. `full` requires
returned bodies and complete source coverage; `partial` retains unverified source
coverage. Missing results/process errors prevent complete publication. A timeout
remains a timeout even after an admission checkpoint. The stopped intermediate
repair runs are excluded.

With OCP 8.0.1.1 and Python 3.12, reproduce the separate reference capture using
its frozen protocol and the same private inputs:

```sh
python3 -B tools/research/occt_census.py collect --sources DURABLE_PRIVATE_SOURCES --selection ORIGINAL_PRIVATE_SELECTION.json --hashes docs/research/step-census/source-hashes.json --out NEW_PRIVATE_REFERENCE --workers 2 --timeout 60
python3 -B tools/research/check_occt_census.py
python3 -B tools/research/compare_import_census.py --native-csv NATIVE_CAPTURE/census.csv --native-manifest NATIVE_CAPTURE/public-summary.json --occt-rows REFERENCE_CAPTURE/rows.json --occt-manifest REFERENCE_CAPTURE/identity.json --out NEW_PUBLIC_COMPARISON
```

The protocol hash is
`8b4a470a00d7afa48481d2cb544e4c7a53797f771ae8cf50717cbb6063918e73`;
its private full reference-row hash is
`815a1bd4320a2a246cc13650cfc2a136fd7f20d67a8ec0b1c5420eb93739eb9b`.
The public reference extract preserves settings profiles and unobserved counts,
and omits source-root/selection paths. The comparison checks every original ID
and source hash, complete terminal coverage, native build identity, category
consistency and the ten retained controls before writing publication artifacts.

The historical admission validation includes 48 scoped release controls, eight observer controls,
six collector/public-path controls, direct OCCT builder checks (including internal healing) on both fixed real workflows,
the Linux format gate and CI-selected eighteen-crate all-target Clippy. Reference
validation includes five actual authored solid/sheet/open-shell/edge/unreadable
controls. Six comparison controls reject changed/missing/duplicate sources,
wrong build attribution, surface-as-solid publication and leaked paths while
preserving partial/timeout unknowns. Public artifacts were independently checked
for counts, hashes and path omission.

The export qualification adds ten integrity controls, two native observer
controls and actual post-builder checks on both authored bodies. Partial manifests,
changed STEP bytes, changed body count and source hash are rejected. An authored
mixed valid/invalid pair cannot qualify the file. The Linux formatting gate and
exact CI-selected nineteen-package all-target Clippy passed. The separate census
audit verified all 3,000 IDs/hashes, the frozen 400/808 population, actual Cargo
build identity, the earlier 322/61/16/1 post-builder totals (including internal healing), paired reference outcomes, ten retained real
controls, unknown timeout work and public path/geometry omission.

The finite-circle repair adds a complete matched-main frozen-400 observation,
32 verified recoveries without losses, five native analytic controls and fifteen
actual STEP/material controls.
Its independent byte audit and per-stage work are retained with the
[matched comparison](step-census-arc-parameters/README.md).

## Sources

- [Current export-validity table](step-census-arc-parameters/results.csv), [manifest](step-census-arc-parameters/results.json), [matched comparison](step-census-arc-parameters/comparison.json), [prior export capture](step-census-export-validity/README.md), [admission table](step-census-periodic/results.csv), [admission manifest](step-census-periodic/results.json), [reference comparison](step-census-periodic/comparison.csv), [reference settings](step-census-periodic/occt-reference.json), [withdrawals from the initial analytic capture](step-census-analytic/initial-analytic-withdrawals.json)
- [Prior 9.20% table](step-census/results.csv), [manifest](step-census/results.json), [source hashes](step-census/source-hashes.json), [retained controls](step-census/retained-controls.json)
- [Original R8 census](../milestones/R8-public-step-census.md), [incidence follow-up](../milestones/R8-incidence-gaps.md), [spline certification follow-up](../milestones/R8-spline-certification.md)
- OCCT [STEPControl_Reader](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V8_0_0/src/DataExchange/TKDESTEP/STEPControl/STEPControl_Reader.hxx) and [ShapeFix_Shape](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V8_0_0/src/ModelingAlgorithms/TKShHealing/ShapeFix/ShapeFix_Shape.hxx)
