# Frozen census through the per-body tier consumer

All **3,000 of the original 3,000 inputs** were observed on 16 Nov 2026 UTC.
The seed is **20261111** and the source hashes, selection and denominator are unchanged.
The [CSV](results.csv) records every original ID, source hash, first refusal, remaining cause, tier, work and export outcome.
The [manifest](results.json) binds the actual producer, Cargo artifact/executable and protocol; it grants this result to the tier consumer, rather than to the earlier source-storage/healer steps.

| File outcome | Count | Share of 3,000 |
| --- | ---: | ---: |
| Complete-source admitted solids | 418 | 13.93% |
| Mixed/partial bodies or incomplete source inventory coverage | 119 | 3.97% |
| Retained source-surface sets, without an admitted B-rep | 2,202 | 73.40% |
| Recognition/geometry-read refusal | 259 | 8.63% |
| Opening timeout | 2 | 0.07% |
| Supported source representations opened, including partial/surface data | **2,739** | **91.30%** |
| Every expected solid exported and valid after OCCT builder processing | **384** | **12.80%** |

The original **186/3,000 (6.20%, about 6%)** was an admission-only result.
The complete-source solid and independent export measures above remain separate from the broader retention count.
An unverified source-surface set is structured source geometry with retained trims and diagnostics; it supplies no exact model, closed material, checked STEP or surface-validity certificate.
The [import contract](../../pieces/08-step-import-repair.md#per-body-opening-tiers-and-source-replay) owns these distinctions and replay semantics.
Existing recipe/application entries do not silently switch to the new consumer.

There are **14,020 recognized source bodies** in opened files: 384 in the valid-solid tier, 663 repaired/approximate solids, and 12,973 surface representations.
Numerical interpretation remains approximate even when retained 3D geometry has not moved.
Actual files were written for 915 bodies: 879 are valid after builder processing and 36 invalid; no reader-refused output body was observed.
Valid exports from a partial file never qualify that entire file.
Export outcomes across all inputs are 384 entirely valid solid files, 104 valid exported subsets, 17 files containing an invalid output, 31 export refusals, 2,202 retained surface sets lacking checked B-rep export, 259 unopened inputs, and three timeouts.
The timeout set includes the two opening timeouts and one export timeout, so it is not added again to the opening-outcome denominator.

## No-loss comparison and actual gains

The baseline is the landed 405 complete solid files: the frozen [400-file admission population](../step-census-periodic/results.json) plus the five files from the [landed cone reference](cone-reference.json).
The export baseline is the previous [369 valid files](../step-census-circular-pcurves/results.json) plus those five valid cone files, for **374**.
The run loses **zero** admitted solid files/body inventories and **zero** export-valid files against those baselines.

Thirteen actual files newly reach complete-source solid admission:
`00000005`, `00000007`, `00000251`, `00000813`, `00001328`, `00001549`, `00001692`, `00001708`, `00001777`, `00002104`, `00002228`, `00004183`, `00004657`.
Ten also export every body valid after builder processing:
`00000005`, `00000007`, `00000251`, `00000813`, `00001549`, `00001692`, `00001708`, `00002104`, `00002228`, `00004183`.
These are observed files, not predicted recoveries from first-refusal counts.
The other three have typed export refusals and receive no export-valid credit.

`00009338` was checked once within this integrated run.
It retains the working numerical whole-file route: one repaired/approximate solid, 1,503,230 opening work including 35,710 source-reading work, and 596,903 aggregate export work.
Its actual output SHA-256 is `473d4b03528bab22b4ca8e20b60bc3864086ab8216cd474632ccc26bc6163df4` and it is valid after builder processing.
No healer attempt is needed on this successful route.
The earlier [healer capture](../tolerant-import-healing/README.md), at its own stated producer, still records the separate atomic face-stage budget refusal; it is not relabelled as this run.

The current cone reference already exports both `00004152` bodies valid at producer `88bb7851a00f0a8531fd85618a2cf083fcda510c`.
The older `63456a28` two-invalid-body diagnosis predates retention of the cylinder V correction and is superseded.
No special rerun, cylinder algorithm change or new source-family credit was taken from that stale diagnosis.

## What remains

Recognition/geometry-read refusals are 66 unsupported torus cases, 54 offset surfaces, 40 solids with voids, 31 linear-extrusion surfaces, 40 operation-budget refusals and 28 malformed/unsupported-root parses.
A refusal means the current reader/query did not produce the source representation under these controls; it is not a finding that the original STEP file is corrupt.
Among retained surface files, remaining cause sets include operation-budget refusal (1,536 files), unsupported numerical families (534), incidence gaps (229), invalid numerical trims (148), unverified intersections (135), and reported self-intersections (85).
These overlapping cause sets are not disjoint file totals.
Curved finite embedding, native-parameter incidence and multi-bound containment also remain explicitly unresolved where the numerical consumer did not establish them.

Thirty export refusals are `step.numerical-gap` and one is `numbers.operation-budget`.
The two opening timeouts are `00000187` and `00008201`; `00004777` opens as a solid but its aggregate export times out.
The fourteen fully admitted files with invalid exports are inherited from the prior frozen qualification; three further invalid-output files have partial results and do not gain full-file credit.
The census does not certify later edits, global self-intersection freedom, material agreement or a complete source-to-result deviation bound.

## Matched inputs, checks and retained parts

Producer: `3437c18d2fcc59bbc74c1a4123bf443237f99f73`.
Actual executable SHA-256: `a9eb07a5ec9e424bcd4463bb4c7c322e845a220a689c5c996317e1122f638f8b`.
Published CSV SHA-256: `39fb11e6bf04113854e2c47dd81582e12e13bce38c52638338ccaee40557d942`.
Selection SHA-256: `1b11410b220d5e8a3156794f6c8f9cc09091a43dcac9af801215fb3a38d571a5`.
Source-hash manifest SHA-256: `3323e104c377023846073a05d7e164ca1082dc3dcf32cc046eee897401454cc7`.
Archive identity: `52e6dd1b6fa38e3cd99af59b662370829129540030975919c3e256dce6ad1dbe`.
Every selected file's actual source digest was checked before and after observation; archive names and private CAD geometry are not published.

A subsequent review fix propagates a caller's deadline through per-body strict construction and checks either outcome before publishing or healing; expiry retains the prepaid source surfaces.
The capture above remains credited to its original producer. Its observer uses deadline-free Deciders and an external process timeout; the follow-up changes no geometry, tolerance, work charges or deadline-free decisions.

One unchanged **20M** allowance covers complete source reading and all opening/admission/healing attempts.
A separate unchanged **20M** allowance covers aggregate export of all admitted bodies.
Observer metadata, serialization, hashes used only for observation, disk I/O and external OCCT work are outside numerical work.
There are two observation processes, one Catallax kernel worker per process, and one 60-second bound per source including its independent export check.
Only after all 3,000 terminal rows exist does the collector publish `complete: true` and percentages.
Incomplete representation coverage retains the frozen source-inspection qualifier at the identical source hash and cannot gain full-source credit.

Independent checks use OCP 8.0.1.1 and the pinned raw-builder helper `5633cdcf16efc13f328ed85b33a738f6a9296e2afb925f4c64bbed7c91025376` for solids.
The corresponding STEP builder overload handles an already admitted surface model.
Builder precision and maximum tolerance are 1e-7 mm, unit factors are millimetres/radians, and the caller disables standard same-parameter and maximum-precision modes as in the prior captures.
The builder may still process/heal internally; `valid after builder` does not mean unhealed validity or source fidelity.
No original CAD file or written artifact is rewritten before that check.

All ten P06 parts retain their original source bytes, source SHA-256, serialized original controls, per-body results and actual outputs in private capture storage.
Their compact outcomes are retained in the manifest.
Validation is limited to formatting, exact hosted scoped all-target Clippy, focused tier/replay/changed-observer tests and this one complete census acceptance run.
The normal 0.001 mm gap cube reopens with 6 faces, 12 edges, 8 vertices, approximate status and identical history, by replay from source plus controls; there is no new healed-body codec.

The existing [OCCT reference](../step-census-periodic/occt-reference.json) uses the same source identities and denominator: 2,993 transferred/opened files and 2,701 healed-valid solid files (90.03%), plus 44 open-shell and 243 other results, one failure and eleven timeouts.
Its default transfer followed by default ShapeFix is a different method and does not establish common tolerance, work or complete source-item coverage.
The comparison is retained as context, rather than rerunning or renaming those reference observations.

To reproduce, build the clean producer with `tools/research/opened_import_census.py build --record <outside-checkout-record>`.
Run its `collect` command with that record, the frozen selected source directory, original private selection, published source-hash manifest and a new outside-checkout capture directory under OCP 8.0.1.1.
The observer header owns its individual-file invocation; `--retain-controls` saves source/control artifacts for the selected ten parts.

## Retained-surface sewing qualification

The single 3,000-file run at kernel producer `b9f7fc355032563d26479215dbd78b0628f9296b` uses two workers, one 20M opening allowance and one aggregate 20M export allowance per file, a 60-second file deadline, and a 6,144 MiB cap on each producer and independent reader.
It records 418 complete-source solid files and 384 files with every solid exported valid, identical file sets to the frozen tier capture above.
There are no new complete solid or export-valid files from this sewing change.
The remaining outcomes are 128 partial files, 2,360 retained-surface files, 92 unopened inputs and two opening timeouts.
The opening count differs from the older producer; this run does not attribute that difference to sewing.

The opening timeouts remain `00000187` and `00008201`; `00004777` retains its completed solid opening and times out only at export.
No observer reaches the memory cap and no process error is recorded.
The collector's stage-attribution correction was checked against the saved rows without rerunning observers.
It reports no measured tier or export loss and no unmeasured baseline pass in the applicable stage.
The original and corrected aggregation receipts remain separate.

The public normal-offset surface fixture also sews into a numerical solid, replays from its immutable source and controls, and writes an actual STEP file.
The pinned OCCT builder reads that artifact as valid with one solid, six faces, twelve edges and eight vertices.
The smaller caller cap retains the source surfaces instead.
Finite skew-line correspondence and complementary-conic-arc controls exercise the added sewing checks.
These controls establish the implemented family, not recovery of all 2,202 older retained-surface files or export authority for arbitrary retained carriers.

## Multiple numerical face-bound qualification

Kernel producer `92a62c83f99afcea4009019b1b6aff4dcdf70d92` was measured once on the same 3,000 source hashes, with the preceding run's two workers, work limits, deadline and 6,144 MiB observer cap.

| Qualified file set | Stitching capture | Multiple-bound candidate |
| --- | ---: | ---: |
| Complete-source solids | 418 | 550 |
| Every solid exported and valid after OCCT builder processing | 384 | 479 |

The measured gains are 132 complete-source solid files and 95 export-valid files.
No previously qualified file loses its measured outcome, and no previous pass is unmeasured in its applicable stage.
All new export-valid files are among the new complete-source solid admissions.
The remaining 37 new admissions gain no export-valid credit.
The first causes of the new solid files were 108 `numerical.invalid-trim`, 16 `self_intersection_not_verified` and eight `not_solid`.
`00005188`, previously blocked by a complete circular outer rim with a hole, now exports valid; `00004869` still remains a surface result.

Other opening outcomes are 129 partial files, 2,227 retained-surface files, 92 unopened inputs and the same two opening timeouts.
Across full and partial files, 43 have `invalid_after_builder` as their primary export outcome and 43 refuse export; 44 contain an invalid written body, including one export-refused file.
The timeout set is unchanged, with `00004777` retaining its completed opening verdict.
No observer reaches its memory cap and no process error is recorded.

The focused controls pass polygonal multiple holes, circular hole pairs, nested/touching/crossing refusals, caller-tolerance clearance and a thin annulus with different rim phases.
Actual double-bore and annulus STEP files also pass the pinned OCCT builder's validity, volume and inside/outside checks.
These remain numerical results with the recorded sampled and continuous-circle checks, not whole-body source-fidelity or exact proofs.

## Numerical planar export qualification

Kernel producer `fd24f13ab8fbc92ba4d5e8eca6ad3d6605b12e90` was measured once on the same 3,000 source hashes, with two workers, unchanged 20M opening and aggregate export allowances, a 60-second file bound and 6,144 MiB caps.

| File measure | Multiple-bound capture | Numerical export candidate |
| --- | ---: | ---: |
| Complete-source solids | 550 | 550 |
| Every solid exported and valid after OCCT builder processing | 479 | 507 |
| Primary `invalid_after_builder` export outcome, full and partial | 43 | 46 |
| Any invalid written body, full and partial | 44 | 46 |
| Export refusal, full and partial | 43 | 1 |

The 28 newly export-valid files were already complete-source admissions, each previously stopped by `step.numerical-gap`.
The 550-file solid set and its body inventories are unchanged; there are no measured losses or unmeasured previous passes.
Of 31 complete-source numerical-gap candidates, three still produce invalid output and receive no credit: `00001687`, `00005501` and `00008723`.
Twelve of the preceding qualification's 13 newly admitted gap cases now export valid; its other 24 invalid-output cases remain invalid.
The difference between primary invalid outcomes and files containing an invalid body accounts for an invalid prefix already written before an export refusal.

The same opening timeouts (`00000187`, `00008201`) and export timeout (`00004777`) remain run facts.
No observer reaches its memory cap and no process error occurs.
The sole remaining export refusal is `numbers.operation-budget` for `00002204`; it is not a geometric diagnosis.

The planar writer retains source vertices, lines and planes and rechecks the selected decimal payload at the body's unchanged tolerance.
It no longer demands exact-zero source incidence, but rejects a reversed finite-line walk and retains existing collapsed-payload refusals.
Fourteen focused STEP/oblique controls pass, including decoding a written source gap without healing it.
Actual `00000835` and `00000018` outputs have all three and all sixteen bodies valid after the pinned builder, respectively.
The first exported `00000835` body also takes a valid material-removing box cut without altering its input bytes.
These observations grant neither a whole-body source-fidelity bound nor export-valid credit without independent checks.

## Cylindrical period-chart qualification

Kernel producer `06910e4076e41530620fb66d8a811c4fdac69297` was measured once on the same 3,000 source hashes, with two workers, 6,144 MiB caps, unchanged 20M allowances and the 60-second file bound.

| File measure | Numerical export capture | Period-chart candidate |
| --- | ---: | ---: |
| Complete-source solids | 550 | 550 |
| Every solid exported and valid after OCCT builder processing | 507 | 516 |
| Any invalid written body, full and partial | 46 | 37 |
| Export refusal, full and partial | 1 | 1 |

All nine newly export-valid files were already fully admitted and had an invalid-after-builder export outcome.
Their IDs are `00000265`, `00000269`, `00001929`, `00004622`, `00004623`, `00005002`, `00005763`, `00007405` and `00008723`.
There are no measured admission/body-inventory losses, export-valid losses or unmeasured previous passes.
`00001687` and `00005501` remain invalid and receive no new credit.
The same two opening timeouts, one export timeout and `00002204` export work refusal remain run facts; no memory termination or process error occurs.

The first invalidity in `00008723` was localised before selecting the fix: the reader moved circular pcurves into a different U period while retaining the seam branches, making one cylindrical face unorientable.
A diagnostic common-period translation restored validity; cyclic loop rotation alone did not.
Independent authored annulus controls reproduced the sensitivity for both face senses.
The candidate passes six periodic STEP tests and ten actual period-gauge/signed-height annulus artifacts with the expected topology and volume.
All three real-part bodies now validate, with unchanged written 3D entity tables.
Only the existing four-sided cylindrical seam-cell chart is translated; 3D geometry, source data, tolerance and the existing numerical checks remain unchanged.
Builder processing can occur internally, and no complete source-fidelity or unhealed-validity claim is made.

## Rational support projection qualification

Kernel producer `7e7ca60efc48830ceb3b540b1a89fd909355418b` was measured once on the same 3,000 source hashes, with two workers, 6,144 MiB caps, unchanged 20M allowances and the 60-second file bound.

| File measure | Period-chart capture | Determinant/trace candidate |
| --- | ---: | ---: |
| Complete-source solids | 550 | 550 |
| Every solid exported and valid after OCCT builder processing | 516 | 516 |
| Any invalid written body, full and partial | 37 | 37 |

There are no measured admission/body-inventory losses, export-valid losses or unmeasured previous passes. The same two opening timeouts, one export timeout and `00002204` export work refusal remain run facts; no memory termination or process error occurs. The candidate receives no additional solid or valid-export credit.

The candidate preserves positive Gershgorin margins and otherwise attempts a positive-definite symmetric projected-Jacobian determinant/trace bound. Thirteen focused tests pass, including a non-ruled cubic sheet without diagonal dominance, indefinite and negative-diagonal refusal controls, and the existing ruled-chart controls. Controls, knots, weights, complete-domain coverage, join checks, tolerances and proof limits stay unchanged. No shared import helper changes.

Six initial refusals change without completing a solid file: `00003966` advances from the support bound to an existing spline endpoint incidence gap; `00005482`, `00006039` and `00007760` advance from that support refusal to a resource limit; `00001166` and `00001707` change from a resource limit to incidence gaps. Resource limits are run facts, not geometric refusals. Five screened real files still retain surface bodies; passing a support proof alone does not admit a solid.

## Capped analytic-family probe (17 Nov 2026)

A separate 12-file probe selects the first four sorted remaining explicit IDs for each of `TOROIDAL_SURFACE`, `OFFSET_SURFACE` and `SURFACE_OF_LINEAR_EXTRUSION`. Its main-fork baseline is `db18482fd9b6d9cb8802cd4f9c633ae05815ca04`; the combined native-tier/family candidate is `1be2ad4028fc803de3572b9f7ce4bb4ecb5fb651`. Identical source hashes, two workers, separate unchanged 20M opening/export allowances, 60-second file bounds and 6,144 MiB process caps are used.

| Capped file measure | Baseline | Candidate |
| --- | ---: | ---: |
| Complete-source solids | 0 | 0 |
| Every solid exported and valid | 0 | 0 |

There are no lost or unmeasured baseline passes. All 12 candidate files retain source surfaces with operation-budget causes. Changed refusals receive no solid or valid-export credit. This probe neither replaces the 3,000-file qualification nor attributes a family-only gain to the new default native tier.

At that candidate, seven focused original-source projection cases pass, including scaled units, a spindle inner sheet and its offset, oblique ellipse extrusion and a signed negative torus offset. A generated normal-offset half-torus reaches the sampled tolerant-solid tier and writes STEP that the direct OCCT 8.0.1.1.0 builder reads valid at 1e-7 mm, without an explicit ShapeFix call. Its native `Geom_OffsetSurface` definition and parameters remain present; volume is 197.39208802180022 mm³ versus 197.39208802178717 expected.

At producer `778bd0ead5c57af305b6654aa1a0ce0f1f0d164c`, the oblique half-ellipse extrusion reaches the same sampled tier and writes an artifact, but its raw-reader UV-period closure is unresolved and receives no export-valid credit. Keeping master-curve flags or rotating the closed walk does not resolve that case. Dropping PCURVE associations also fails raw validity and changes the observed volume, so it is not an accepted workaround. No complete-source body-deviation or unhealed-validity claim is made.
