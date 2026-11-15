# B03 bounded inventory: required fallback remains unmapped

## Families needed beyond B02

No additional family is demonstrated necessary for a named Phase 1 production
request. This is an evidence blocker, not a no-consumer closure or port completion.

| Family | Actual caller and shipped fixture | Evidence / required status |
| --- | --- | --- |
| Zero-result tracer fallback on OFF | `DropTracer::trace_curve2` in `surface_drop_curve.rs`; `tests/fixtures/drop-curve/cases.txt` physical line 64, `drefuse`; generator names `OFF` | Executed here: zero general-drop curves, reached fallback, typed refusal. Phase 1 production-request necessity unproved. |
| Failed or multiple-result general-drop fallback | Same helper; no specifically mapped required fixture identified | Source-supported branches only; no new execution or port justification. |
| Projected side search | `drop_and_trim_curve(..., on_surface=false)`; physical line 65, `BUMP.iso(0,1.3)` | Executed typed refusal; separate unsupported branch, not tracer fallback. B01 scope decision remains. |
| General closed seam / parameter remapping | Seam refusal in `surface_drop_curve.rs`; B02 guards and helper witnesses in `surface_trace.rs` | B02 covers one bounded native-parameter family. Remaining refusals do not justify another family or a remap; unchanged here. |

The named OFF fixture is identified but cannot be mapped to a required Phase 1
consumer from the inspected public evidence. Per the assignment, publish this
bounded inventory instead of implementing speculative cases. Its complete
positive pinned donor result is also not shipped: the saved gate stores `i-1`
and its README describes the actual answer as privately retained. No private/raw
output was accessed or regenerated.

## Baseline and method

Actual upstream main was observed, fetched and frozen on 2026-11-15 UTC to
`48a6d53d57300d7b248e45a3aa86ebca6ad7f06d`, the merged B02 PR #607 commit.
Fresh branch: `phase1/B03-general-drop`. No prior local/remote B03 branch or B03
implementation PR was found. Prior P04/B01/B02/N02 work and evidence are preserved.
No applicable AGENTS.md or .agents/skills/SKILL.md was present. The checkout
verification skill was read for context; no CLI qualification workflow invoked.

Read the [P04 map](../P04/PORT-CONSUMERS.md), merged B02 guards, known-miss list,
public drop fixture README/generator, all 14 [P02 questions](../../../evidence/phase1/P02/m10.jsonl),
their comparison and both published worker-1/worker-4 run/operation records.
Static JSON extraction does not execute recipes. New geometry execution was
Rust only, on the unchanged shipped drop inputs.

Production architectural route: `pair/trace_candidates.rs` ->
`Intersector::do_intersection` -> side-range drop in `surface_trace.rs` ->
`drop_and_trim_curve`. Curve/surface and line/surface coincidence also reach it
through `coincident_ranges`. This establishes possible reachability, not that a
particular required request enters the missing fallback. The already-ported
`surface_drop_general.rs` walker is distinct from the missing
`SmSurfaceTracer::TraceCurve` / `FlushCurve` fallback. Candidates never replace
independent Catallax acceptance checks.

## The 14 recorded public P02 requests

Historical P02 outcomes, not fresh current-main results. All retain limit 20M.
The operation logs fail at STEP for 001–003, 007 and 008, and Boolean for the
others. None records internal drop-branch telemetry. Export refusal or budget
exhaustion cannot be reclassified as reached fallback; an unsupported-pair
refusal does not establish this branch was reached or that all earlier branches
were untouched.

| ID | Recipe | Historical refusal | Used work |
| --- | --- | --- | ---: |
| M10-001 | side-subtract | `step.unsupported-geometry` | 8,219,519 |
| M10-002 | side-union | `step.unsupported-geometry` | 8,219,568 |
| M10-003 | side-intersect | `step.unsupported-geometry` | 8,214,023 |
| M10-004 | top-subtract | `numbers.operation-budget` | 19,999,916 |
| M10-005 | top-union | `numbers.operation-budget` | 19,999,916 |
| M10-006 | top-intersect | `numbers.operation-budget` | 19,999,916 |
| M10-007 | through-cut | `step.unsupported-geometry` | 2,527,978 |
| M10-008 | m7-cut | `numbers.geometric-gap-unproven` | 1,802,159 |
| M10-009 | import-subtract | `intersect.parametric-coincidence` | 115,984 |
| M10-010 | import-union | `intersect.parametric-coincidence` | 115,984 |
| M10-011 | import-intersect | `intersect.parametric-coincidence` | 115,984 |
| M10-012 | spline-pair-subtract | `intersect.unsupported-pair` | 74,693 |
| M10-013 | spline-pair-union | `intersect.unsupported-pair` | 74,693 |
| M10-014 | spline-pair-intersect | `intersect.unsupported-pair` | 74,693 |

[Structured inventory](../../../evidence/phase1/B03/p02-inventory.json) retains
all complete final refusal/work records and original refusing-operation records
at both published worker counts. It references complete question/response/operation
files with hashes. Recipe names are resolved through the portable wrappers,
including the import wrappers; no tally is used as necessity evidence.
No P02 request was rerun and no existing result/expected answer changed.

## Executed fixture observation

[Focused Rust test](../../../crates/catallax-trace/tests/b03_drop_inventory.rs):
only the two unchanged shipped `drefuse` rows. OFF is the cubic with controls
(0,0,.35), (1,.5,.55), (2,-.5,1.05), (3,1,.95), on the generator's tilted plane,
interval [0,1], tolerance .001, keep-all true, on-surface true. The other row
requests on-surface false.

A temporary Rust diagnostic observed `DropTracer::trace_curve2 curve_count=0`
for OFF. Direct general drop with keep-all false also returns zero curves,
max-drop 0 and max-deviation 0. Those zero fields on an empty candidate result
are not geometric zero-distance evidence. The wrapper returns the complete
`NotPorted(Intersect)` refusal. The projected row exits the separate false-branch
guard with the same refusal. Diagnostic production bytes were restored in a
`finally` block and checked equal to baseline; the committed test reran without
instrumentation. No production algorithm is changed or translated; no new donor
adaptation or positive differential result is claimed.

The focused test passes both unchanged integer refusal expectations and compares
no positive binary64 result. See [validation](../../../evidence/phase1/B03/VALIDATION.md).
It is not a finite-curve completion or new geometric answer.

## Exact blockers and reopening requirements

- A named authorized Phase 1 request, or bounded extracted public input with
  production-caller provenance, dynamically reaching the required fallback beyond
  B02; distinguish failure, zero curves or multiple curves and retain domains.
- Complete public pinned differential outputs for that required family: finite
  curves/UVs, statuses and deviations, with known-miss controls. The current
  refusal sentinel cannot qualify a positive port. Oracle work belongs to its
  separately authorized lane.
- If a remap is needed, establish it only for that named consumer with interior
  3D/UV correspondence; otherwise preserve B02's typed refusal.

B03 remains blocked. This does not prove the fallback unnecessary or promise it
will resolve the 14 refusals. No new curve acceptance, proof, STEP writer,
E01/N03 work, tolerance/cap, saved answer, hosted workflow, C++, crash/sanitizer
probe or private input is involved. No new application-work or worker-count
qualification is claimed. Earlier captures remain historical evidence.
