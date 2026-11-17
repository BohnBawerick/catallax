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
| Smooth-roof/top-box proxy, not the frozen housing | 20,220,641 | `fillet.numerical-body`: one support is planar |
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
The proxy's later edit was not attempted after its fillet refusal. The actual
frozen curved-housing recipe is `tools/checkpoint/phase2/curved-housing.json`.
The loft regression preserves the bounded-refusal contract. Successful edits
must keep face, edge and vertex names and pass independent geometry checks.
See the [numerical method and final-admission split](research/spline-fillet-offset-search.md#16-november-numerical-multi-patch-body-and-raw-step).

Historical cases retained for phase 5: the 1e-4 mm carrier/cap-only request
took 29,335,696 work, and exact offset proof across atlas seams on
perpendicular planes exhausts at 100M. Optional numerical-body handover
encoding also exhausts 100M; no transfer-engine performance change was made.

## Part-two radius candidate, 17 November

The original 0.1 to 0.15 mm radius change was localized to the first end-plane
cut at 100M. A bounded-interpolation candidate kept the residual threshold,
feature tolerance and work allowance unchanged, but did not finish the edit.

| Stage | Charged work | Result |
|---|---:|---|
| Multi-patch loft fillet body | 92,314,860 | Built; independent radius/cap geometry passes |
| Independent initial-body radius/cap checks | 626,608 | Passed in a separate proof scope |
| Radius edit, 1/10 to 3/20 mm | 99,999,973 | `numbers.operation-budget` in the first end-plane cut |

The candidate passes the existing two-chart carrier/cap control and all three
admission controls, but still does not build the radius edit. No edited-name,
edited-geometry or edited-STEP success is claimed. The candidate was parked
and the original bisection retained. The listed rows are separate measured
request scopes, not an internal profile or one combined request. Internal
charges before the failing cut were not collected. Steering on 17 November
defers another radius attempt and moves on to import and the real frozen
housing.

## Imported-part continuation, 17 November

The retained 27,210-byte unfilleted loft STEP is read through the public import
API, before numerical embedding and filleting. Its SHA-256 is
`859a5f0f1443cf394ece5bb6d7c07dcb79d6de5ff79bbc36560110f03303411a`.

| Stage | Charged work | Result |
|---|---:|---|
| Import application | 98,427,317 | `numbers.operation-budget` |
| Numerical embedding, admission and fillet | Not reached | No imported fillet result |

All charged work belongs to `import.apply`. The retained cause is
`SelfIntersectionNotVerified` for source solid `#21`, faces `#23` and `#207`,
because its verification exhausted the work allowance. This is not a proved
self-intersection. No import implementation, tolerance or work limit changed.

## Frozen housing and M8 continuation, 17 November

The actual housing is
`tools/checkpoint/phase2/curved-housing.json`, SHA-256
`803558dbb48dab4add6fc1ba70d3ea8a2ac11ba044b8256d2c7ea9030c51e7e6`.
Its recorded edit changes the middle section radius from 14 to 15. Each
chain uses one 100M request. Incremental checkpoints retain only prefixes
already charged in that same request.

| Variant | Loft work | Tool work | Cut work | Total work | Result |
|---|---:|---:|---:|---:|---|
| Original | 1,371,690 | 50 | 139,526 | 1,511,266 | `intersect.unsupported-pair` |
| Middle-radius edit | 1,371,698 | 50 | 139,526 | 1,511,274 | `intersect.unsupported-pair` |

Both stop at `cut`: rational atlas with Bezier extrusion spline, outside
piece-4 pair coverage. Neither reaches the fillet, so neither establishes
an editable filleted housing or its STEP. These are the real arch-housing
results, not the earlier smooth-roof/top-box proxy.

All 13 original files under
`crates/catallax-features/tests/fixtures/p5-later-cuts/` were then run
unchanged, each under its own 100M request:

| Case | Prefix work | Final step work | Total work | Result |
|---|---:|---:|---:|---|
| k01 wall hole | 1,058,398 | 98,941,443 | 99,999,841 | `numbers.operation-budget` |
| k02 top slot | 1,058,398 | 1,224,882 | 2,283,280 | `intersect.parametric-surface`, non-separated transverse crossing |
| k03 rim notch | 1,058,398 | 24,344,499 | 25,402,897 | `boolean.parametric-query`, finite solver unavailable |
| k04 torus pocket | 1,058,398 | 98,941,599 | 99,999,997 | `numbers.operation-budget` |
| k05 plane notch | 1,058,398 | 82,013,695 | 83,072,093 | `numbers.operation-budget` |
| k06 cylinder drill | 1,058,398 | 82,010,092 | 83,068,490 | `numbers.operation-budget` |
| k07 cone drill | 1,058,398 | 98,923,535 | 99,981,933 | `numbers.operation-budget` |
| k08 sphere dimple | 1,058,398 | 98,931,318 | 99,989,716 | `numbers.operation-budget` |
| k09 ring groove | 1,058,398 | 81,575,808 | 82,634,206 | `numbers.operation-budget` |
| k10 axis section | 1,058,398 | 81,573,472 | 82,631,870 | `numbers.operation-budget` |
| k11 offset section | 1,058,398 | 82,009,922 | 83,068,320 | `numbers.operation-budget` |
| k12 far disjoint | 1,058,398 | 9,038,752 | 10,097,150 | Built, 22 faces |
| k13 fillet first | 21,973 | 3,449,241 | 3,471,214 | `intersect.unsupported-pair`, cone/torus sign unproved |

For k13 the failing final step is `drilled`; the other rows reach `later-cut`.
Only k12 builds, and its separate STEP request writes in 260,103 work,
SHA-256 `203ae123b49186c571eb47fb5f0b8b5859c3e95f8bbc2214f27e073154d6c891`.
A disjoint subtraction adds no intersecting-cut coverage. All stages over
20M remain phase-5 cases; no work limit or tolerance was widened.
