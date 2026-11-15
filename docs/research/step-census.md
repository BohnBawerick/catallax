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

The complete repaired 15 November kernel census imports **378 of the same 3,000
files with full source coverage: 12.60%, beside 9.20% at the prior checkpoint and
the original 6.20% (about 6%)**. All 3,000 terminal observations retain the original
seed, selection, individual source-byte hashes, denominator and observation caps.

| Outcome | 11 November | Prior 15 November | Repaired analytic import |
|---|---:|---:|---:|
| Imported with every source item covered | 186 (6.20%) | 276 (9.20%) | 378 (12.60%) |
| Imported, some source items not verified | 3 (0.10%) | 3 (0.10%) | 3 (0.10%) |
| Refused with a named cause | 2,782 (92.73%) | 2,719 (90.63%) | 2,617 (87.23%) |
| External observation timeout | 29 (0.97%) | 2 (0.07%) | 2 (0.07%) |
| Total observed | 3,000 | 3,000 | 3,000 |

All 276 prior full imports and all three partial imports remain in those
categories. **102 additional files import fully**, each previously stopped at
`planar_repair_not_verified`. The current categories are 190 certified imports
(including the three partial imports) and 191 numerical approximate imports.
Compared with 11 November, 191 previously refused files and one previously timed
out file now import fully. The manifests list the actual model IDs.

The [initial analytic capture](step-census-analytic/initial-results.json)
[reported 382 full imports](step-census-analytic/initial-results.csv) (12.73%). After restoring
independent source-pair checks, four of those files now return
`self_intersecting`: `00001233`, `00005142`, `00007233`, `00007433`.
[Their identities and causes](step-census-analytic/initial-analytic-withdrawals.json)
are retained; the initial 12.73% is a historical checkpoint, not the current
result. That historical capture records its clean collector revision and supplied
binary digest; it predates the actual Cargo-artifact binding now required for the
current run. Two intermediate repair runs stopped before completion and contribute
no complete census claim.

The [current per-model table](step-census-analytic/results.csv) retains each
original first refusal and the prior 9.20% outcome/cause beside the new outcome.
The [current manifest](step-census-analytic/results.json) binds the clean producer,
actual Cargo artifact, observer, collector, toolchain, caps and table hash.
The [prior table](step-census/results.csv) and [manifest](step-census/results.json)
remain unchanged.

### First causes before and after

| Returned cause | 11 November | Prior 9.20% | Repaired 12.60% |
|---|---:|---:|---:|
| Planar source repair not verified | 0 | 975 | 129 |
| Source incidence gap | 0 | 505 | 502 |
| Unsupported source type or variant | 1,314 | 876 | 846 |
| Resource limit | 334 | 126 | 171 |
| `not_solid` | 203 | 80 | 80 |
| Self-intersection not verified | 85 | 79 | 78 |
| `self_intersecting` | 822 | 53 | 101 |
| Malformed input | 24 | 24 | 24 |
| Numerical family unsupported | 0 | 1 | 645 |
| Numerical finite trim rejected | 0 | 0 | 41 |
| Total returned refusals | 2,782 | 2,719 | 2,617 |

These are first stops, not predicted admission gains. Every exact unsupported
variant and transition remains in the manifests. Initial source bytes are now
charged before kernel verification hashing, including mismatches. This exposes
some resource limits earlier; a smaller unsupported bucket does not establish
new STEP-type support. Actual admission gain over 9.20% is 102 files.

The finite analytic fallback retains shared 3D vertices/carriers and reconstructs
finite parameters for planes, cylinders, lines and circular arcs. Cheap source
checks run across every solid before optional curved proofs. The latter share
one collective 1M-work prepayment inside the unchanged 20M request. Completed
checks and every unresolved pair remain explicit; neither unresolved coverage
nor zero sample residuals grants a global proof. Numerical STEP remains marked
approximate. Full periodic bands and other unsupported families still refuse.
Two follow-ups remain: with less than 1M work remaining a file can refuse at the
collective prepayment, and the legacy planar intersection predicate is not fully
metered. Recorded counters cover existing reservations and added traversal/setup
checkpoints, rather than every legacy planar operation.

All 191 numerical imports use the fixed 0.01 mm tolerance. Their endpoint and
curve/face sample checks do not establish continuous between-sample separation
or a complete source/body bound. Existing local source-to-derivative bounds and
new pointwise residuals retain separate columns; an absent source bound is not
zero. This census measures kernel admission. The
[analytic qualification](../phase1/import-analytic-numerical/README.md) checks
actual import, pose and raw STEP output on two fixed real files; the
[earlier oblique qualification](../phase1/import-oblique-step/README.md) covers
its own fixed planar sources. No edit/export claim is made for every admission.

### Standard OCCT comparison on all 3,000 files

The complete reference run uses `STEPControl_Reader` normal root transfer,
`ShapeFix_Shape::Perform` with default settings, then `BRepCheck_Analyzer`, via
OCP 8.0.1.1. Source bytes are checked before and after every observation. The
[frozen protocol](../../tools/research/occt_census.py),
[reference manifest/default profiles](step-census-analytic/occt-reference.json),
[per-file comparison](step-census-analytic/comparison.csv), and
[comparison manifest](step-census-analytic/comparison.json) retain its scope.

| Exclusive reference outcome | Files | Share |
|---|---:|---:|
| Valid whole healed result, containing solids and no free lower-dimensional topology | 2,701 | 90.03% |
| Open-shell result with no solid | 44 | 1.47% |
| Other non-null opened result | 243 | 8.10% |
| Failed | 1 | 0.03% |
| External timeout | 11 | 0.37% |
| Total observed | 3,000 | 100% |

Opening a non-null transferred shape was observed for **2,993 files**, including
five that later timed out. An opening checkpoint does not establish completed
healing, validity, complete source coverage or a solid. The other six timeouts
have no observed opening. Root transfer counts and whole-shape/solid checks are
retained separately. Empty cells mean unobserved, not a successful check.

All 378 Catallax full-import source files produce reference valid-solid results.
Another 2,321 files with a Catallax refusal and both Catallax timeouts produce
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
Both current runs read the durable private source corpus. Every source hash is
checked before and after its observation. The
[ten fixed P06 capability sources](step-census/retained-controls.json) remain
retained privately with original bytes and hashes; all ten still import fully.
Only public model IDs, hashes, diagnostic scopes and aggregate metadata are
published, under the [existing licence scope](../milestones/R8-public-step-census.md#source-and-licence-scope).

The current native producer is
[`d78cb705d85a6d98391a5469a4e086cf353804e8`](https://github.com/BohnBawerick/catallax/tree/d78cb705d85a6d98391a5469a4e086cf353804e8).
Its actual Cargo observer artifact SHA-256 is
`28bf7cc68d757e4333193dba1ab4187520460deba977de07565309a64dc819b0`.
The current per-file CSV has SHA-256
`d6f6d6674a70d4a7506b73859e5d875164b6047e5b420839308d09c9d476d7a5`.
The source/body fallback is numerical and opt-in at the kernel API. SolidDraft's
numerical-group and replacement-reference UI work remains with the app slice.

From that clean pinned producer, run:

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

Native validation includes 48 scoped release controls, eight observer controls,
six collector/public-path controls, raw OCCT checks on both fixed real workflows,
the Linux format gate and CI-selected eighteen-crate all-target Clippy. Reference
validation includes five actual authored solid/sheet/open-shell/edge/unreadable
controls. Six comparison controls reject changed/missing/duplicate sources,
wrong build attribution, surface-as-solid publication and leaked paths while
preserving partial/timeout unknowns. Public artifacts were independently checked
for counts, hashes and path omission.

## Sources

- [Current native table](step-census-analytic/results.csv), [manifest](step-census-analytic/results.json), [reference comparison](step-census-analytic/comparison.csv), [reference settings](step-census-analytic/occt-reference.json), [withdrawals from the initial analytic capture](step-census-analytic/initial-analytic-withdrawals.json)
- [Prior 9.20% table](step-census/results.csv), [manifest](step-census/results.json), [source hashes](step-census/source-hashes.json), [retained controls](step-census/retained-controls.json)
- [Original R8 census](../milestones/R8-public-step-census.md), [incidence follow-up](../milestones/R8-incidence-gaps.md), [spline certification follow-up](../milestones/R8-spline-certification.md)
- OCCT [STEPControl_Reader](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V8_0_0/src/DataExchange/TKDESTEP/STEPControl/STEPControl_Reader.hxx) and [ShapeFix_Shape](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V8_0_0/src/ModelingAlgorithms/TKShHealing/ShapeFix/ShapeFix_Shape.hxx)
