# Phase-5 work cases

The application allows 100M until phase 5; the target remains 20M. These completed
edits exceed the target and retain their measured stage work for later speed work.
No optimisation is part of their capability qualification.

| Case | Result | Work used | Evidence |
|---|---|---:|---|
| M10-004 | builds at 100M, over the 20M target | 26,777,693 | [Stage work and raw STEP checks](phase1/spline-boolean-top-cuts/REPORT.md) |
| M10-005 | builds at 100M, over the 20M target | 27,289,073 | [Stage work and raw STEP checks](phase1/spline-boolean-top-cuts/REPORT.md) |
| M10-006 | builds at 100M, over the 20M target | 25,945,726 | [Stage work and raw STEP checks](phase1/spline-boolean-top-cuts/REPORT.md) |
| 2.3 spline sheet, trim, extension, +2/5 mm thickening and STEP | builds at 100M, over the 20M target | 29,145,161 | [Signed STEP and healing-off OCCT receipts](SURFACE-BODIES.md#reproduction-and-observed-proof) |
| 2.3 same chain, -2/5 mm thickening | builds at 100M, over the 20M target | 29,142,981 | [Signed STEP and healing-off OCCT receipts](SURFACE-BODIES.md#reproduction-and-observed-proof) |

The 2.3 counts cover each complete construction/export request; individual stage
counts were not profiled. Both signed solids pass raw OpenCASCADE validity with
healing off. Profiling belongs to phase 5.

The slow loft-source admission test timed out at 1,200 seconds in CI. It stays
available with `--ignored`. The complete 100M top requests also have an explicit
slow regression. Required Python CI still checks all three 20M refusals and artifact
absence, in the `step-regressions-slow` job. The listed full-request counts come from the recorded source before the
review's narrower predispatch checks; those long requests were not repeated.

## Incomplete 2.1 bounded attempts

Each row is a separate request capped at 100M, not evidence that the combined
body-build/export chain fits one 100M app request. The charged counters below
are reported as emitted; `exhausted` requests can report fewer than 100M units.
All rows above 20M remain phase-5 work cases.

| Stage | Charged work | Result |
|---|---:|---|
| Multi-patch loft fillet body | 96,732,996 | Built |
| Direct numerical STEP and reread | 80,265,942 | Built; raw OpenCASCADE read and later cut pass |
| Radius edit, 1/10 to 3/20 mm | 99,999,999 | `numbers.operation-budget` |
| Imported-part fillet request | 98,434,067 | `numbers.operation-budget` during input import/evaluation |
| Housing loft, spline cut, fillet | 20,220,641 | `fillet.numerical-body`: one support is planar |
| M8 k01 wall hole | 99,999,841 | `numbers.operation-budget` |
| M8 k02 top slot | 2,283,280 | `intersect.parametric-surface`: non-separated transverse crossing |
| M8 k03 rim notch | 25,402,897 | `boolean.parametric-query`: finite curve solver unavailable |
| M8 k04 torus pocket | 99,999,997 | `numbers.operation-budget` |
| M8 k05 plane notch | 83,072,093 | `numbers.operation-budget` |
| M8 k06 cylinder drill | 83,068,490 | `numbers.operation-budget` |
| M8 k07 cone drill | 99,981,933 | `numbers.operation-budget` |
| M8 k08 sphere dimple | 99,989,716 | `numbers.operation-budget` |
| M8 k09 ring groove | 82,634,206 | `numbers.operation-budget` |
| M8 k10 axis section | 82,631,870 | `numbers.operation-budget` |
| M8 k11 offset section | 83,068,320 | `numbers.operation-budget` |
| M8 k12 far disjoint | 10,097,150 | Built, 22 faces |
| M8 k13 fillet first | 3,471,214 | `intersect.unsupported-pair`: cone/torus sign unproved |

Milestone 2.1 remains open: stable radius-edit names, an imported-part fillet,
the completed housing chain and intersecting later M8 cuts are not proved.
The housing's later edit was not attempted after its fillet refusal.
The permanent loft regression checks geometry and either stable edited names
or a typed work refusal; it does not claim a successful radius edit.
See the [numerical method and final-admission split](research/spline-fillet-offset-search.md#16-november-numerical-multi-patch-body-and-raw-step).

Historical cases retained for phase 5: the 1e-4 mm carrier/cap-only request
took 29,335,696 work, and exact offset proof across atlas seams on
perpendicular planes exhausts at 100M. Optional numerical-body handover
encoding also exhausts 100M; no transfer-engine performance change was made.
