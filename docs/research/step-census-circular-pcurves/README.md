# Constant-height circular STEP pcurves

The complete frozen-population comparison records **369 of the original 3,000
files (12.30%) valid after OCCT builder construction, which can heal pcurves
internally**. The matched integrated control has 361 valid files. **Eight actual
files recover; every previous valid file remains valid.** All ten retained P06
controls remain valid. The original seed, source hashes and denominator remain.

| Frozen 400, after builder processing including possible internal healing | Matched control | Constant-height repair |
| --- | ---: | ---: |
| Every expected body exported and valid | 361 (90.25%) | 369 (92.25%) |
| An exported body invalid | 22 | 14 |
| Export refused | 16 | 15 |
| Timeout | 1 | 2 |

The [paired comparison](comparison.json) records the eight recovered files:
`00000718`, `00002553`, `00003124`, `00005223`, `00005497`, `00005774`,
`00009338`, `00009656`. All were previously invalid, and there are **zero losses**.
File `00002204` changes from aggregate export work refusal to export timeout;
it earns no recovery credit. The other timeout is `00004777`.

The [3,000-ID table](results.csv) retains first refusals and after-fix causes.
All 400 imports reproduce the expected 808 bodies. The other 2,600 inputs were
not reevaluated for new admissions. Original admission-only results, including
186 (6.20%, about 6%), remain in [the census](../step-census.md).
The separately retained OCCT result is 2,701 of 3,000 (90.03%) valid after normal
source transfer and default ShapeFix healing, with different precision settings.
Neither measure establishes unhealed validity, full source fidelity or a
common-tolerance comparison.

## Reproduction and bounded family

Real part `00009338` reproduced an open finite cylindrical wire. A written
circular pcurve had a tiny varying height, inherited from source vertices that
were slightly off the circle plane. During builder processing that map shifted
by a full period relative to neighbouring edges. The source body and its
canonical shared-vertex trim walk had already passed numerical admission.

The [original authored baseline](authored-incidence-baseline.json), producer
`326ceeed02613d05e42acade5258c957a8424fff` with the prior STEP kernel, reproduced
144 invalid outputs among 224 native-admitted cases. It covers signed vertex
height errors from `1e-14` to `1e-10` mm, four equivalent charts, four vertices,
and cardinal/oblique poses. Zero-error controls remain in the matrix.

The repair follows OCCT's [circle projection on a cylinder](https://github.com/Open-Cascade-SAS/OCCT/blob/V8_0_0/src/ModelingData/TKGeomBase/ProjLib/ProjLib_Cylinder.cxx).
A circle whose normal is parallel to the cylinder axis has constant axial
height. The writer uses the height of its written centre when a finite circle's
existing pcurve has varying height and the written rational directions pass a
cheap parallelism check. It retains the finite U chart, all existing flat maps,
tilted carriers and the canonical shared-vertex walk used for admission.
Independent checks on the actual rounded pcurve still sample nine positions,
including both endpoints, against the retained source tolerance.

Only actual pcurve coordinates change. All written 3D geometry, vertices,
surfaces and oriented loops remain. Numerical STEP remains approximate, with no
complete body bound. The [3D graph audit](three-dimensional-graph-audit.json)
compares all 699 declared written body files, including unchanged pcurve surface
associations and preference. It allows no loop rotation or reversal. Thirty-three
byte streams change, and every compared 3D graph and source-tolerance label agrees.

Before implementation I inspected OCCT's curve projection, STEP edge/loop
translation and `XSAlgo_ShapeProcessor::CheckPCurve`, together with the
[P04 SMLib map](../../phase1/P04/PORT-CONSUMERS.md). The remaining SMLib seam/tracer
omissions concern intersection work. This repair changes only the Plane/Cylinder
numerical STEP exporter; cone export, geometry reading, spline intersection,
Boolean paths and shared topology helpers are unchanged.

## Regression controls and measurement identity

An earlier [broad endpoint-projection experiment](withheld-broad-projection.json)
completed all 400 observations at 363 valid files, with four gains and two losses
(`00001784`, `00005348`). It was withheld. It introduced tiny slopes into flat
maps of tilted carriers. Both losses had previously flat affected pcurves.
The [32 tilted-carrier controls](authored-withheld-tilt-baseline.json), producer
`7a5bb9c909ac4ba33dca09051cbae673244d2f55`, reproduced eight invalid outputs under
that broad experiment. The restricted repair preserves both lost real files
byte for byte against the integrated control.

The measured candidate producer/kernel is
`2d53b8b1bc307b075c6ff8323aed0595a9e93468`, executable SHA-256
`8b7e4f17a7ec3c0cedc3ba940d14c6fb9708a2b4de6e105d06399f45cc57212e`.
Its exact [CSV](results.csv) has SHA-256
`bbe34caeb3161804261583637de29d5ef4597f9c1918bcb510b56d643c1dd142`;
the [manifest](results.json) binds its actual Cargo artifact and protocol.
Later publication/rebase heads are distinct from this measured producer.

The matched control producer/kernel is
`1e8b24b651c468e8f81cc1748fcc6a8b3fe3db3b`, executable
`2d917b760606d9ba8d3aa94ff502936aaed5b51320b247430c59c4af8af49962`.
The [control manifest](matched-control.json) and [audit](matched-control-audit.json)
record its complete run after the intervening main integration. Its CSV is
byte-identical to the [previous 361-file table](../step-census-cylinder-seams/results.csv),
SHA-256 `890808462b91320ae054d226462ad34b04bf378dbced33a2e81fc57795e67974`.
All 699 STEP files also reproduce the previous capture. Candidate and control
import work, source hashes, categories and body counts agree for all 400 files.
Their only changed production path is `crates/catallax-step/src/numerical/analytic.rs`.

The unchanged [protocol](../step-census-export-validity/README.md#what-qualifies-a-file)
explains OCCT's internal healing and the legacy `raw_*` names and false
observer-healing flags in preserved schemas. Those caller-level settings do not
establish unhealed validity, and no per-file internal-healing count is available.

## Validation and limits

Work-PC validation spans 15–16 November 2026 UTC, Rust 1.99.1, OCP 8.0.1.1.
Twenty-two native tests pass, including source vertex/curve tilt over-tolerance,
rounding, distant-chart, collapsed-arc and exhausted-work refusals.
All [256 quarter-cylinder controls](authored-controls.json) pass validity,
expected volume/topology and 768 material probes. The previous 224 full-cylinder
STEP controls/448 probes and 15 finite-circle controls/45 probes also pass:
495 authored STEP files and 1,261 material probes in total.

Generate the new controls outside the checkout with
`CATALLAX_ANALYTIC_EVIDENCE=PRIVATE/authored cargo test --release -p catallax-step --test numerical_cylinder_incidence`.
Replay them in the documented OCP environment with
`python tools/research/check_circular_pcurves.py PRIVATE/authored --out PRIVATE/qualification.json`.
The [focused real controls](focused-real-controls.json) check both former losses
byte for byte and the named recovery before collecting another census.
The [source-reference control for 00009338](retained-real-9338.json) retains its
source hash, one source root and one returned solid; both volumes are
32,799.79710011322 cubic millimetres and all 45 material probes agree.

The [independent audit](independent-byte-audit.json) verifies the actual Cargo
artifact/executable, all 3,000 IDs/hashes and the denominator, 400 durable source
hashes, 808 expected body counts, every declared STEP hash and millimetre/radian
units. Of 699 declared body files, 677 reach builder checks: 645 valid bodies and
32 invalid bodies. A valid body inside an invalid file earns no file credit;
the 808 imported bodies are not a claim that all exported successfully.

| Stage | Observed counters | Maximum work | Sum of observed work |
| --- | ---: | ---: | ---: |
| Import | 400 | 18,929,210 | 571,316,874 |
| Aggregate export of each file's bodies | 398 | 14,150,281 | 202,430,646 |

Each stage retains its separate unchanged 20M allowance. Both export timeouts
have unknown export work. Their missing counters prevent a performance inference
from the smaller observed maximum or sum. The whole-file deadline remains 60
seconds, with two file workers and one kernel worker each. OCCT and observer
I/O/hashing remain outside kernel counters. All 15 returned refusals are
`step.numerical-gap`. Timeouts and remaining failures stay in the published table.

Owner-requested [certified-wire triage](certified-wire-triage.json) inspected
`00002036` from admission through direct source/export builder diagnostics.
Its 16 source wires are single positive-radius closed circles with exact endpoint
incidence. Source and export share the direct-builder periodic-band failures;
normal source transfer adds periodic topology and matches the independently
integrated rounded-tube volume. No 3D source-wire crossing was demonstrated for
this file; the other eight certified failures remain unqualified. This triage
reused the frozen capture and did not recollect the population.

These checks do not qualify arbitrary edits, spline admission or SolidDraft's
replacement UI. The [cause observations](cause-observations.json) retain the
remaining invalid files without predicting further recoveries.
