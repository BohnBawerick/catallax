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

Independent OCP 8.0.1.1 checks read both the actual imported source and the
emitted k12 STEP through raw `StepToTopoDS_Builder`, with healing,
SameParameter and TransferRoots off. Both are valid closed single solids
and take a valid, material-removing reference cut:

| File | Faces | Volume before, mm³ | Volume after, mm³ |
|---|---:|---:|---:|
| Unfilleted source STEP | 10 | 8.057291666666668 | 6.516691814464737 |
| k12 far-disjoint STEP | 22 | 12558.682045617654 | 11544.082343087084 |

The restored permanent loft regression passes with the original bisection:
96,732,996 work for body admission, 626,466 for independent geometry, and a
typed radius-edit refusal at 99,999,999 in its separate request. No successful
strict radius edit is inferred from the passing bounded-refusal contract.

## Native rolling-ball continuation, 17 November

The new standalone native cylinder-rim continuation exhausts an unchanged
100M request at 1e-7 mm, or 0.1 nm, both before and after separating binary64 analytic
proposals from interval acceptance. The refusal is
`numbers.operation-budget` inside `NumericalFilletSpine::walk`, before its
independent geometry controls. The failing request's exact final work count
was not retained. This is an incomplete new continuation path, not a work
diagnosis for the owner's everyday fillet refusals.

The same run builds the cone rim at 13,724,976 work and checks the
three-plane common ball and singular-contact refusal at 5,746. Successful
cone continuation retains 259 checked stations. No budget or tolerance was
widened. Further walker work follows the ordered trimmed-clearance repairs.

## Four-contact fan material queries, 17 November

The R5 concurrent four-plane apex builds in 230,473 work and writes STEP in
1,607,549. Its four B-rep containment requests exceed the 20M target but stay
within the unchanged 100M application allowance:

| Point, mm | Material | B-rep query work |
|---|---|---:|
| (0, 0, 7.9) | Inside | 28,360,368 |
| (0, 0, 8) | Outside | 27,267,790 |
| (2, 0, 7) | Inside | 53,190,147 |
| (2, 0, 7.8) | Outside | 35,775,093 |

Field compilation uses 14,474 work; each corresponding field query uses
257,785. An earlier test incorrectly shared one allowance across compilation
and all eight independent queries. All four controls now pass with a separate
unchanged 100M allowance per public request; none was dropped or weakened.
Raw OCP 8.0.1.1 independently agrees with all four controls at 1e-7 mm and
checks a valid closed single solid and a material-removing edit, with healing,
SameParameter and TransferRoots off. Query cost remains a phase-5 case.

## Ordinary numerical fillet application request, 17 November

`catallax_app::build_incremental_without_display` now retains the cross-chart
loft fillet as a typed numerical outcome. Building its exact prefix and
resuming the fillet share one unchanged 100M allowance: 96,733,000 work.
The declared tolerance remains 1e-7 mm. Canonical face origins and an
independent fitted-geometry control pass; no exact checkpoint is produced.

The existing direct numerical proof used the release profile. The first
application check instead used an unoptimized test build and exceeded its
900-second SSH deadline; its final work and outcome were not recovered.
The comparable release application regression passes in 169.69 seconds
(305.90 seconds including compilation). No kernel algorithm, tolerance or
work cap was changed to obtain this result. Future remote proofs retain
launch PIDs, logs and a wall deadline covering the entire offload group.

After the main rebase and app-owned outcome cutover, the same whole
prefix-plus-fillet request passes at 96,731,167 work in 143.72 seconds,
still above 20M under the unchanged 100M allowance and 1e-7 mm tolerance.
This is not a simple native fillet or an internal stage profile. The
separate repaired-body/export-provenance consumer also passes.

After the generic feature-result payload, whole-curve UV retargeting and
verified embedding integration, this ordinary consumer passes at 96,863,834
whole-request work in 148.57 seconds (255.01 with compilation). It still
uses 100M and 1e-7 mm. The preceding bounded attempt timed out after
1056.73 seconds in the test, before its owner diagnostic ran; the approved
exact disjoint-UV-box predicate is included in the passing integration.

## General rolling-ball contact stage, 17 November

The new material-side/radius regression constructs all 12 untrimmed contact
strips of a 20 × 30 × 40 mm box. Each row is one complete request with an
unchanged 100M allowance and 1e-7 mm tolerance, including guides, contact
walks, fitted strips/rails and independent geometric controls:

| Radius, mm | Request work |
|---|---:|
| 0.5 | 52,541,470 |
| 1 | 26,868,969 |
| 2 | 15,604,461 |
| 5 | 9,139,129 |

R0.5 and R1 exceed the 20M target and remain phase-5 cases. These are not
closed-body fillets: neighbour-domain trimming, corner assembly, ordinary
application wiring and owner acceptance are not included. No internal
stage profile was collected, so the radius-dependent cost has no attributed
cause. The release test passed in 59.09 seconds; its bounded command,
revision, process-group/PID records, full log and exit are retained in
`proof-general-contact-box-2/` in the task evidence directory.

The separate nonplanar-neighbour request passes in 33,629,407 work at the
same allowance and tolerance. Its scope includes one R1 strip, construction
of the quadratic neighbour `z = 10 + x²/100 + y²/50`, a selected transverse
cut and independent world/both-pcurve geometry controls. It is not an
isolated profile of the cut solver. Actual trimmed neighbour-face domains
and closed-body admission are not included. The release test takes 1.11
seconds; `proof-general-neighbour-cut-1/` retains the bounded command,
revision, PID/group records, full output and exit.

## Capped rolling-strip solid, 17 November

One R1 edge of a 20×30×40 box now completes the general contact walk,
two actual neighbour intersections, retained-border/support-loop edits
and unchanged numerical solid admission at 1e-7 mm. The whole request,
including its original prefix and numerical source conversion, spends
72,898,257 work under the unchanged 100M allowance. It is above the 20M
target and remains recorded for phase 5; no internal stage profile was
collected. The independent six-plane area controls run after that request.

The release test takes 27.54 seconds (100.91 with compilation). This proves
one closed capped-strip body, not network/corner publication, STEP exchange
or the owner part. No cap, tolerance or admission guard was changed.

## Shared rolling-strip and corner network, 18 November

The corner owner's boundary-helper proof includes all 24 paired section
curves, rotated sphere charts and independent controls. These whole requests
are not component profiles or admitted bodies:

| Radius, mm | Request work |
|---|---:|
| 0.5 | 80,909,382 |
| 1 | 63,408,116 |
| 2 | 57,956,045 |
| 5 | 56,609,707 |

The integrated R1 request separately publishes all 12 strips and eight
corners of a 20×30×40 box as one admitted numerical solid. It spends
72,340,370 work, including the source, contacts, paired boundaries and full
body admission. Independent area controls check all six original plane
domains after that request. The release test takes 159.83 seconds, or
240.21 seconds with compilation. These costs remain above the 20M target
under the unchanged 100M allowance and 1e-7 mm tolerance. No stage profile
or optimisation was added. Owner, application and STEP acceptance are not
included in this body proof.

The first ordinary six-edge R25 owner-rim request at R1 reaches general
rolling construction and refuses `shapes.rational-atlas-unproven` while
admitting the fitted strip's regular chart cover. Whole prefix-plus-fillet
work is 75,128,530, not exhausted. The command takes 175.49 seconds including
compilation. This is a failed owner request, not a passing fillet or an
isolated chart-proof cost.

## Final native source embedding admission

The embedding owner's final P5 source proof at
`f16d9d680c8b81feffc86df05610d5f24901332b` admits all 41 original edges
of its 22-face body, retaining 11 added paired seam edges. Its unchanged
tolerance is 1e-7 mm and its whole operation allowance is 100M.

| Stage | Work |
|---|---:|
| Source construction | 1,073,248 |
| Native embedding | 21,468,886 |
| Native B-rep admission | 29,943,064 |
| Whole request | 52,485,198 |

Embedding, admission and the whole request remain over the 20M target.
These are the owner's measured stage counters, not a CPU profile or a
like-for-like timing benchmark. This proves source admission, not a fillet,
radius edit, later cut or STEP exchange. Owner, demo, bracket and P1 source
admission also pass below 20M in that handoff.

## Owner radius overflow eligibility

At `fc1205299ef397eac024cd199eef39cdc303ebea`, the retained R25 owner
source, six-edge selection and R5 rolling request exhaust 99,999,988 work
under the unchanged 100M allowance and 1e-7 mm tolerance. Runtime is 44.94s
(51.63s command). This is a failed whole request, not an isolated stage
profile or a passing fillet.

The initial cheap contact proof accepted only Line carriers. A separate
source-only diagnostic spends 516,540 work and identifies the thin-wall
edge as an Intersection curve on native cylinder and parallel-plane
supports. The geometric constraint therefore needs no Line representation.
At `3b138022e22f599169b9e1ce2a7912295ce45158`, the corrected constraint
uses the native supports rather than the edge's representation:

| Request scope | Work | Result |
|---|---:|---|
| Source, six-edge selection, early R5 contact proof | 535,459 | R4 and the actual input side face; 0.05s |
| Ordinary application source/checkpoint and six-edge R5 call | 670,379 | Typed R4 refusal, input face #57, terminal hint |
| Source and one curved owner strip at R3.99 | 61,539,101 | Numerical strip and independent radius controls; 255.80s |

The first two requests are below 20M; the ordinary application command takes
59.12s including 54.39s standalone compilation. The R3.99 cost remains above
the target. It is one strip, not the six-edge rim or an admitted fillet body;
the full-rim 100M failure remains open. Walk and fit components were not
separated in this measurement. No tolerance, physical guard or cap changed.

At `b4a3dc39676d698ff944e7eeeb734e6764185435`, noncharging test
snapshots separate that same 61,539,101-work R3.99 request:

| Measured component | Work | Retained stations | Completed station checks |
|---|---:|---:|---:|
| Offset-centre walk | 17,888,503 | 324 | 7,686 |
| Carrier and rail fit | 42,872,488 | 324 | 8,655 cumulative |

Other source, selection and guide setup was not separately measured.
The smallest walked tau span is approximately
[0.9989062718, 0.9993796174], width 4.733456089e-4.
The fit component remains above 20M and is assigned to its owner; the
walk and actual split-cylinder walk refusal remain separate correctness
investigations. These are request work snapshots, not a CPU profile.
No production geometry, accounting, guard, tolerance or allowance changed.
