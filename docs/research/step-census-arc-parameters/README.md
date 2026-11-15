# Finite-circle STEP parameter comparison

The complete comparison records **354 of the original 3,000 files (11.8%) valid
after OCCT builder construction, which can heal pcurves internally**, beside
the retained **2,701 of 3,000 (90.0%) valid after normal OCCT transfer and default
ShapeFix healing**.
Neither method establishes unhealed validity, source fidelity or a
common-tolerance comparison.

The frozen population remains 400 admitted files and 808 expected bodies.
The original seed, all source hashes and the 3,000-file denominator are unchanged.
The other 2,600 inputs were not reevaluated for new admissions.
Original admission-only counts, including 186 (6.20%, about 6%), remain in
[the census](../step-census.md).

| Frozen 400, after builder processing including possible internal healing | Matched main | Finite-circle repair |
| --- | ---: | ---: |
| Every expected body exported and valid | 322 (80.5%) | 354 (88.5%) |
| An exported body invalid | 61 | 29 |
| Export refused | 16 | 16 |
| Timeout | 1 | 1 |

The [paired comparison](comparison.json) identifies **32 actually recovered
files, with every previous valid file retained and no regressions**.
All ten retained P06 controls remain valid after builder processing.
The [cause observations](cause-observations.json) distinguish the 44 files with
observed whole-period UV jumps from the 32 verified recoveries; a diagnostic
group size is not a predicted gain.
Remaining failures stay in the result table.

## Isolated change

A STEP `EDGE_CURVE` carries no native parameter range.
OCCT reconstructs its range by projecting the vertices.
The old writer could build pcurves using a different periodic branch, leaving
finite wires open in UV although their 3D endpoints coincided.
Real part `00000292` reproduced this failure before the change.

The analytic numerical writer now anchors each finite circle's natural lower
endpoint at phase zero before constructing its pcurves.
It retains the written center, axis, radius, endpoint identities, finite arc
direction and face UV charts; only the circle reference phase changes.
The actual rounded frame must establish a positive reference ray, at most
`1e-12` radians of phase error, and at most the smaller of `1e-9` mm and the
configured scalar rounding limit in phase displacement.
An unresolved phase or finite range is refused.
Existing source tolerances and numerical edge/surface checks still apply.
Numerical results and STEP remain approximate, with no complete body bound.

This follows OCCT's treatment of omitted ranges and projected periodic
parameters in [STEP edge translation](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_0/src/DataExchange/TKDESTEP/StepToTopoDS/StepToTopoDS_TranslateEdge.cxx)
and [periodic parameter adjustment](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_0/src/DataExchange/TKDESTEP/StepToTopoDS/StepToTopoDS_GeometricTool.cxx).
Anchoring avoids a crossing of the carrier's periodic origin rather than
relying on a reader to choose the same branch.
The [P04 SMLib map](../../phase1/P04/PORT-CONSUMERS.md) places its remaining
seam-start and tracer work in curve/surface intersection.
This exchange repair needs neither tracer and changes no spline intersection,
Boolean, native admission or shared numerical topology helper.

## Measurement identity

The candidate producer is `3651cd367cbb2f5d28118b06c6431193eb1b0d62`.
Its actual Cargo executable SHA-256 is
`17b928c6ee2bb18e522f3676a7db52fdfd64f2bb6dc6b323f137c87a3b2917bc`.
Its exact [result CSV](results.csv) has SHA-256
`0f9c1727baa480e4c4c9e446e9495fcfa891dfff0249a7acfbddbcc3e6928431`.
The [capture manifest](results.json) retains build and protocol identities.
Later publication and rebase heads are distinct from this measured producer.

The matched control producer is `a207ea676f52ee077cbb186674407dea3e64e1d0`, with
production kernel `f0abed8ed0ec5b8cc66626d60141bc69e99e5823`.
Its executable is `1c304cdcc090508bd8d9a533c37ffd6009263a0d862c7529e64e4b1019b4446f`;
its [CSV](matched-main.csv) is
`8e765971a02100100e4c2a1bb7b81046bd9b2d01d763f7a1ab8b13baf318be3f`.
The [control manifest](matched-main.json) records the actual build.
The only production difference between this control and the candidate is
`crates/catallax-step/src/numerical/analytic.rs`.

Main contained other kernel changes since the older 322-file producer, so a
matched control was necessary before attributing recoveries.
The control reproduced every earlier outcome and emitted STEP hash; only its
import work counters differed from the older capture.
Candidate and matched-control import work, source hashes, categories and body
counts agree for all 400 files.
The source protocol and both captures remain byte for byte.
The [existing protocol](../step-census-export-validity/README.md#what-qualifies-a-file)
explains OCCT's internal healing and the legacy `raw_*` names and false
observer-healing flags in these unchanged capture schemas.
Those fields do not establish unhealed validity.

## Validation and limits

Validation on 15 November 2026 UTC used Rust 1.99.1 and OCP 8.0.1.1.
Five native analytic tests pass, including twelve equivalent periodic charts,
pose rebuild, unchanged distant-chart and exhausted-work refusals, and a
four-digit oblique phase refusal with a passing seventeen-digit control.
The [authored qualification](authored-controls.json) checks fifteen actual STEP
files, their topology and volume, and 45 material probes.
Fresh and frozen fixture bytes agree.
The other numerical STEP, periodic and oblique suites also pass.

Generate authored fixtures outside the checkout with
`CATALLAX_ANALYTIC_EVIDENCE=PRIVATE/authored cargo test --release -p catallax-step --test numerical_analytic`.
Then run `python tools/research/check_arc_parameters.py PRIVATE/authored --out PRIVATE/qualification.json`
in the documented OCP environment.
Full frozen-population collection uses the existing protocol and clean producer
build records, with separate output directories for control and candidate.

The [independent audit](independent-byte-audit.json) checks the 3,000 IDs/hashes,
400 durable source hashes, all 808 expected body counts, actual Cargo artifact
binding, every declared STEP hash and millimetre/radian units.
It checks 699 written body files; 677 reached builder validation, with 613 valid
bodies and 64 invalid bodies.
A valid body inside an invalid multi-solid file never earns file credit.
The 808 imported bodies are not a claim that all 808 exported successfully.

| Stage | Observed counters | Maximum work | Sum of observed work |
| --- | ---: | ---: | ---: |
| Import | 400 | 18,929,210 | 571,316,874 |
| Aggregate export of each file's bodies | 399 | 19,585,351 | 221,994,048 |

Each stage retains its separate unchanged 20M allowance.
The whole-file native-plus-builder observation retains its 60-second deadline
and two file workers with one kernel worker each.
OCCT work and observer I/O/hashing remain outside the kernel counters.
The 16 refusals remain 15 `step.numerical-gap` and one aggregate export work limit.
File `00004777` again times out during export after import completed; its export
work and resulting validity remain unknown, not zero.
All 400 observations terminated; no partial run is labelled complete.

Builder validity and these authored probes do not prove a complete source/body
deviation bound, material agreement, an arbitrary edit or a Boolean for every
census input.
The SolidDraft replacement-reference UI remains with the app slice.
