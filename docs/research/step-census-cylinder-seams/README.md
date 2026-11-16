# Cylindrical STEP seam ordering

The complete frozen-population comparison records **361 of the original 3,000
files (12.03%) valid after OCCT builder construction, which can heal pcurves
internally**. The matched control has 354 valid files. **Seven actual files
recover, with no losses.** All ten retained P06 real controls remain valid.

The seed, source hashes, denominator and frozen 400 admissions are unchanged.
All 400 imports reproduce the expected 808 bodies. The other 2,600 inputs were
not reevaluated for new admissions. Original admission-only results, including
186 (6.20%, about 6%), remain in [the census](../step-census.md).

| Frozen 400, after builder processing including possible internal healing | Matched control | Seam ordering |
| --- | ---: | ---: |
| Every expected body exported and valid | 354 (88.50%) | 361 (90.25%) |
| An exported body invalid | 29 | 22 |
| Export refused | 16 | 16 |
| Timeout | 1 | 1 |

The [paired comparison](comparison.json) retains the actual recovered IDs:
`00000363`, `00001857`, `00004333`, `00005897`, `00008682`, `00009319`,
`00009321`. They belong to the observed wire-closure group, which contained 19
files before this fix. That group size was not treated as a predicted gain.
The [cause observations](cause-observations.json) retain the remaining failures.

## Reproduction and isolated repair

Real part `00000363` imports with two expected bodies. Its outer cylinder
exported an invalid UV wire with the seam's two pcurves assigned to opposite
uses. A transient OCCT slot-swap diagnostic restored validity and volume;
that modified shape earned no census credit.

The authored reproduction then separated chart phase from loop order:
[28 baseline cases](authored-baseline.json) have both surface-axis signs and
seven equivalent chart phases for a solid cylinder and square stock with a
through bore. All 14 bore cases passed; all 14 outer-cylinder cases failed.
Rotating the closed outer-cylinder loop restored its finite material, while
exchanging the serialized pcurve list or changing preferred representation did
not. An independently generated OCCT cylinder also passed the same builder.

OCCT's [STEP loop reader](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_0/src/DataExchange/TKDESTEP/StepToTopoDS/StepToTopoDS_TranslateEdgeLoop.cxx)
installs both seam pcurves at each occurrence, using cumulative edge, bound and
face orientation. The final occurrence can therefore determine the stored pair.
The exporter now chooses a cyclic start whose final occurrence of each retained
cylindrical seam has forward cumulative orientation. STEP edge `same_sense`,
including descending natural parameters, participates in that calculation.
A walk with no such cyclic start keeps its existing admitted ordering.

Only the start of the closed walk changes. The edge sequence is never reversed;
vertices, carrier geometry, parameter values, pcurve lists and trim directions
remain intact. This follows the reader's explicit orientation rule and OCCT's
[ordered-wire export](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_0/src/DataExchange/TKDESTEP/TopoDSToStep/TopoDSToStep_MakeStepWire.cxx).
The [P04 SMLib map](../../phase1/P04/PORT-CONSUMERS.md) places the remaining
seam-start/tracer omissions in intersection work. No new projection, tracer,
intersection, Boolean, cone reader or shared topology helper is needed here.
The only production change is the analytic numerical STEP writer.

All existing source tolerances, decimal checks and admission refusals remain.
Numerical bodies and STEP remain approximate, with no complete body bound.
The [exchange-graph audit](exchange-graph-audit.json) independently compares all
699 declared written body files: 21 byte streams change, but their headers and
entity graphs are identical after allowing only cyclic `EDGE_LOOP` starts.
It preserves oriented-edge senses and the order of pcurve lists.

## Measurement identity

The measured candidate producer/kernel is
`7563b40c48229cfd57407351f716c3a01867d2a2`. Its executable SHA-256 is
`58128174b490ab39d1aa2927251a4a3769f90dcc5f91528befbb815497b8be25`.
The exact [3,000-ID result table](results.csv) has SHA-256
`890808462b91320ae054d226462ad34b04bf378dbced33a2e81fc57795e67974`.
The [capture manifest](results.json) retains the actual Cargo build identity.
Publication and later rebase heads are distinct from this measured producer.

The matched control producer is `0902b99f3e30893ce11cfb3d69488d55f93b2d92`,
with production kernel `621cb33775218515db99a6399065b4abe8261e54` and executable
`e1f02dac8fdc56636ee3de4166d3a16e1af9d09b15a2ca79327190a454f7b4c1`.
Its [manifest](matched-control.json) records a fresh build and complete run.
Its CSV is byte-identical to the retained [354-file result table](../step-census-arc-parameters/results.csv),
SHA-256 `0f9c1727baa480e4c4c9e446e9495fcfa891dfff0249a7acfbddbcc3e6928431`.
This control incorporates the intervening main integration and typed finite-arc
refusal fix. Candidate and control import work, source hashes, categories and
body counts agree for all 400 files.

The unchanged [protocol](../step-census-export-validity/README.md#what-qualifies-a-file)
explains internal OCCT healing and the legacy `raw_*` names and false observer
healing flags in the preserved capture schemas. They describe caller-level
settings and do not establish unhealed validity. No per-file internal-healing
intervention count is available.

## Validation and limits

Work-PC validation on 15 November 2026 UTC uses Rust 1.99.1 and OCP 8.0.1.1.
Eighteen native tests pass across the cylinder, analytic, periodic and oblique
suites, including the existing rounding, distant-chart and exhausted-work
refusals. The [224 authored files](authored-controls.json) cover both material
sides, both surface-axis signs, seven phases, both seam parameter directions and
all four cyclic starts. Each passes whole-solid validity after builder
construction, expected topology and volume, and two material probes (448 total).
The existing 15 finite-circle STEP controls and 45 material probes also pass.

Generate the cylinder files outside the checkout with
`CATALLAX_ANALYTIC_EVIDENCE=PRIVATE/authored cargo test --release -p catallax-step --test numerical_cylinder_chart`.
Replay them in the documented OCP environment with
`python tools/research/check_cylinder_seam_order.py PRIVATE/authored --out PRIVATE/qualification.json`.

The [actual two-solid control](retained-real-cylinder.json) preserves the source
hash of `00000363`, both source roots and both returned solids. Their volumes
match the normal source-transfer reference within `1e-6` cubic millimetres:
7,853.981633974483 and 7,771.018366025515. All 90 deterministic material probes
agree. This is a finite control, not a complete source-fidelity bound.

The [independent byte audit](independent-byte-audit.json) binds the actual Cargo
artifact and executable, checks all 3,000 IDs/hashes and the denominator,
rechecks all 400 durable source hashes and 808 expected body counts, and checks
every declared STEP hash and millimetre/radian units. Of the 699 written body
files, 677 reach builder checks: 625 valid bodies and 52 invalid bodies.
A valid body inside an invalid multi-solid file earns no file credit.
The 808 imported bodies are not a claim that all exported successfully.

| Stage | Observed counters | Maximum work | Sum of observed work |
| --- | ---: | ---: | ---: |
| Import | 400 | 18,929,210 | 571,316,874 |
| Aggregate export of each file's bodies | 399 | 19,586,631 | 222,004,344 |

Each stage retains its separate unchanged 20M allowance. The whole-file
observation retains its 60-second deadline and two file workers with one kernel
worker each. OCCT work and observer I/O/hashing remain outside kernel counters.
The 16 refusals remain 15 `step.numerical-gap` and one aggregate export work
limit. `00004777` again times out during export after import completes;
its export work and validity remain unknown.

The separately retained OCCT reference is 2,701 of 3,000 (90.03%) valid after
normal source transfer and default ShapeFix healing. It uses different
precision/healing settings. Neither measure establishes unhealed validity,
source fidelity or a common-tolerance comparison. This frozen run does not
qualify arbitrary edits, spline admission or the SolidDraft replacement UI.
