# P04 — remaining port calls and Phase 1 consumers

## Decision supported by this report

Steps 1–9 are landed in the frozen main ancestry, with tracked Rust units and saved differential gates. That means the historical dependency-order slices landed, not that every SMLib path works or any Phase 1 milestone is qualified. The three named omissions are still explicit `NotPorted(Intersect)` returns.

Projected drop has no demonstrated current production caller: all tracked production `drop_and_trim_curve` calls pass `on_surface=true`. Seam starts and the tracer fallback are reachable conditions inside that existing true branch, used for surface-intersection range answers. They are relevant risks to Phase 1.3 and indirectly 1.5, but public evidence inspected here does not attribute any of the 14 Boolean requests, loft top cut or housing cuts to those exact conditions. Do not schedule three implementations solely from an old port line count.

The B01–B03 card definitions are not present in the tracked sources searched at this baseline. If those cards mean projected drop, seam start and fallback respectively, the provisional disposition is: projected drop needs a named consumer before implementation; seam start and fallback need bounded consumer-trigger evidence and missing oracle outputs before implementation. This is a dependency recommendation, not an invented card assignment or a cancellation decision. Firstmate must supply the actual card definitions and choose their priority.

## Baseline, scope and method

- Repository: `BohnBawerick/catallax` only. Remote `main` observed with `git ls-remote origin refs/heads/main`, then fetched and frozen to `f4db345f2e7c51835d1ebeba9f6446d117144117` on 2026-11-14 UTC. Later remote movement does not change this report's baseline.
- Isolated worktree `/workspace/p04-catallax`, branch `phase1/P04-port-consumers`, created from that exact commit. Original `work` checkout and older branches preserved.
- Read tracked sources, fixture READMEs, saved text rows and generator descriptions; inspected first-parent commit ancestry for landed claims. No donor source repository, companion repository, private corpus, retained raw oracle directory, C++ build/execution, hosted gate or CI reads were used.
- No `AGENTS.md` or `.agents/skills/SKILL.md` was found in the repository; `.agents` contained no files. The verification skill in `.claude/skills/verify-catallax/SKILL.md` was read for context; this report does not perform CLI/STEP qualification and does not invoke its full qualification workflow.
- Source hashes and exact saved-fixture row counts are in `evidence/phase1/P04/source-manifest.json`. Rust saved-gate execution is recorded separately in `evidence/phase1/P04/VALIDATION.md`. Historical measurements quoted below remain source reports, not new oracle observations.

All source paths below are relative to the repository at the frozen commit. Line references identify the baseline, not an API contract.

## Call / fixture table

| Remaining path | Exact call boundary and trigger | Current callers and Phase 1 consequence | Existing saved evidence that suffices | Evidence still missing before a positive port claim |
| --- | --- | --- | --- | --- |
| Projected drop | `surface_drop_curve.rs:1262–1273`, `drop_and_trim_curve(..., on_surface=false)` returns immediately. File header names projected side search. | All production direct callers use true: `surface_trace.rs:2459` (side-range drop onto each surface), `curve_surface_global.rs:1226` (coincident-range classification), and `coincident_ranges` at `surface_drop_curve.rs:1329` (used by `line_surface.rs:736`). Only the public test input currently requests false. No present Phase 1.3 blocker is established. A future approximate projection/import/fillet consumer would need its own contract; Phase 1.1/1.2 do not authorize guessing one. | `drop-curve/cases.txt` physical line 65 encodes the bump iso-curve with `on_surface=0`, expected `i-1`. This suffices to preserve refusal. `curve-on-surface` gate covers `EvaluatePoint`, UV/normal evaluation and the global solver's projected branch predicate; it does not cover the missing end-to-end search. | Public complete pinned-oracle output for the false row and bounded positive/empty/multiple-piece projected-search cases; actual consumer input and requested domains. The current row deliberately stores refusal rather than the donor answer. Additional derivative/cache/lifting operations cannot be inferred complete from the predicate gate. |
| Closed-surface seam start | `DropTracer::find_boundary_start_points`, `surface_drop_curve.rs:1007–1013`: exactly two start solutions plus closure in either axis => refusal. Not every closed surface or every seam crossing enters this branch. | Reachable through the true range-drop calls above, so potential Phase 1.3 periodic spline/cylinder/torus side-range blocker; Phase 1.5 inherits intersection needs. No public per-consumer trace inspected proves a named target reaches it. | `surface-drop` point-drop fixtures include two seam answers; `drop-curve` direct `general` cases include cylinder helices across/from seams, seam rulings and u/v orientation; these suffice for their existing operations only. `drop-start` covers 41 singular-start point-value cases, not this closed-seam branch. | A bounded `DropAndTrimCurve(on_surface=true)` fixture that demonstrably hits the two-start refusal, in both closure directions and restricted domains, with complete stored UV pieces and side ownership from an authorized local oracle. No dedicated positive saved gate for this branch was identified. General `DropCurve` seam cases cannot stand in for `DropTracer` seam-start selection. |
| General tracer fallback | `DropTracer::trace_curve2`, `surface_drop_curve.rs:1204–1222`: failed general `drop_curve` maps to refusal, or number of returned curves differs from one => refusal. The source labels missing `SmSurfaceTracer::TraceCurve` / `FlushCurve`. Zero curves also enters the latter condition. | Same true range-drop route; long coincident sides become intersector range curves. This can stop branch candidates for Phase 1.3 and consequently fillet work in 1.5. This is not the already-ported general `SmSurface::DropCurve`. | `drop-curve/cases.txt` physical line 64, tilted cubic OFF, expected `i-1`, preserves the bounded refusal. `general` rows exercise the successful general walker. `intersector` range rows and `curve-surface` coincident ranges cover current successful integrations. | Complete public pinned-oracle output for the refused OFF request; separate safe fixtures demonstrating failure versus zero/multiple-piece fallback and final flushed UV curves, errors and intervals. Existing saved refusal cannot compare the omitted successful donor result. No local oracle was run here. |

The two `drefuse` rows are generated in `tests/fixtures/drop-curve/generate.py:140–149`: OFF is a cubic 0.1 off the tilted plane at one control point; the second is `BUMP.iso(0, 1.3)` with `on_surface=0`. `tests/drop_curve.rs:68–90` accepts `NotPorted` only for refusal operations; a success on those rows fails the gate. Thus “132 cases pass” includes two preserved refusals and is not 132 completed geometric drops. The raw positive SMLib answers are described as privately retained and are unavailable within this task's allowed evidence.

## Phase 1 consumer map

| Consumer | Visible route/evidence | What the three remaining calls establish |
| --- | --- | --- |
| 1.3 spline Boolean cuts, loft top cut, later housing cuts | `ROADMAP.md` 1.3; `docs/M8-LATER-CUTS.md` forced-trace boundary; `catallax-intersect/src/pair/traced.rs:358` -> `Candidates::generate` -> `trace_candidates.rs:134` `Intersector::do_intersection`; fitted plane proposals use `pair/traced/fitted.rs:13` -> `Candidates::plane_section`. Long side ranges in `surface_trace.rs` use the drop helper on each surface. | Seam/fallback are source-reachable risks. Projected false branch is not requested by these calls. Existing per-case refusals include proof budgets, fitted-boundary queries and unsupported pairs; port completion cannot be claimed to remove all of them. No new Boolean replay or STEP export was executed. |
| 1.5 fillets between spline faces across patches | `ROADMAP.md` says it depends on 1.3; `M27-SPLINE-FILLET-SCOPE.md` separates offset-search proof from crossing intersections; `surface_drop_general.rs` already contains the general drop, while `rolling_ball.rs` / `fillet_intersector.rs` are separate candidate units. | Dependence is architectural/source inference. Neither 41 singular-start cases nor general-drop cylinder walks qualifies multipatch fillets, radius-edit naming, topology or export. A fillet consumer trace is needed to assign any of these omissions as its immediate blocker. |
| 1.1 numerical result contract and 1.2 real-file import | `ROADMAP.md` 1.1/1.2; N01 owns the result API under this delegation. | No direct false-branch caller or named import case was established by the inspected call search. Result/approximation policy must be settled by its owner; P04 changes no Result adapter or tolerance. |

Candidate tracing can fail before `prove_selection` / `prove_sections`. A successful trace also still needs the independent complete component cover and exact selector checks in `trace_candidates.rs`. Its finite binary64 conversion and native proposal windows are not proof of completeness, permitted approximation, or accepted geometry.

## Landed steps 1–9: bounded evidence

Every commit below appears in the frozen main's first-parent history. The listed source and saved gate files are tracked at this baseline. “Landed” refers to those scoped units; excluded trimmed-face, directed, fallback and singular cases remain excluded where documented. Gate letters retain the M25 case-driven meaning.

| M25 step | Landed ancestry evidence (abbreviated SHA / PR) | Tracked source and saved gate evidence |
| --- | --- | --- |
| 1 foundation / NURBS helpers | `dd2fc8515e18` / #388 | `scalar.rs`, `vector.rs`, `extent.rs`, `geometry.rs`; `tests/fixtures/foundation/{README.md,cases.txt}` |
| 2 general search / solutions | `c2e0267f2da2` / #396 | `global_solver.rs`, `solution.rs`; `search` fixture directory |
| 3 spline objects | `d2e599fd0180` / #397 | `spline.rs`; `spline` fixtures (README reports 1,075 cases, maximum 3 ULP; not Gate A/B completion) |
| 4 curve services / cache / solve | `d1e46b265be8` / #416 and `818e3cd89a98` / #421 | `matrix.rs`, `local_solve.rs`, `curve_cache.rs`, `curve_solve.rs`, `curve_global.rs`; `local-solve`, `curve-cache`, `curve-solve`, `curve-global` fixtures |
| 5 surface evaluation / cache / point solve | `412dae089341` / #422 | `surface_eval.rs`, `surface_cache.rs`, `surface_solve.rs`, `surface_global.rs`; `surface-eval`, `surface-cache`, `surface-solve` fixtures; no general trimmed-face claim |
| 6 curve/surface, Gate A | `71a76df3e175` / #423 | `curve_surface.rs`, `curve_surface_global.rs`; `curve-surface` fixtures; later range follow-ups do not imply every directed operation |
| 7 marching dependencies | `fa27dc48948b` / #426, `8175c24c5bf9` / #427, `4b0975c86a49` / #428 | `plane.rs`, `surface_drop.rs`, `surface_surface.rs`; `plane`, `surface-drop`, `surface-surface` fixtures |
| 8 intersector, Gate B | `75065d365bdb` / #432, `0d794025ef20` / #433; output prerequisite `73c0f5c20d87` / #429 | `surface_intersector.rs`, `surface_trace.rs`, `surface_curves.rs`; `intersector` and `output-curves` fixtures |
| 9 curved singularity / spline pair, Gate C cases | `ff71ef7e9cd7` / #434 | `surface_trace.rs` Dupin/singularity handling; curved/curved point and intersection cases in `intersector` fixtures; not blanket four-parameter proof acceptance |

Follow-ups include general curve drop `8119912` / #470 and analytic/boundary recovery `6675005` / #579, as well as #498 analytic lines, #508 circular arcs and #563 polar boxes recorded by the known-miss document. These are existing work, not proposed B01–B03 implementations. `docs/THIRD-PARTY.md` is the attribution inventory; some earlier rows retain historical refusal descriptions, so current code plus the known-miss status takes precedence over treating every older row as current.

## Saved gates that can be reused, and their limits

Existing saved public inputs suffice for Rust regression checks of the landed units: foundation, search, spline, local solve, curve services/cache, surface evaluation/cache/point solve, curve/surface, surface/surface, general drop, analytic curves and the intersector. They require no C++ execution. Keep the existing candidate comparison allowances exactly as written in the test harnesses. A partial branch gate is reusable for that branch; it is not an end-to-end answer gate.

For the three omissions, the existing suites suffice to preserve refusals and check neighboring behavior, but not to turn a refusal row into a positive donor comparison. `drop-start` is a private-tracer Rust unit gate under `surface_drop_curve.rs`, distinct from the public integration `surface-drop` point gate. `curve-on-surface` explicitly lacks its own cache; generated branch inputs carry donor node data. It cannot qualify a newly assembled projected-search tree.

The `drop-curve` README records unexercised general-drop seam-jump correction, pole dog-leg straightening, end-of-domain rescue, rescue of a kept far drop and filtering several curves when far drops are not kept. Its surviving mutations include the tangent-angle threshold and boundary-tangent correction. These are limitations even on the already-ported walker. A new local oracle request must target safe defined behavior and preserve full outputs; this task neither executes nor requests memory-fault/sanitizer probes. If such a probe is thought necessary, report that need to Firstmate separately.

## Known misses remain independent acceptance obligations

`docs/M9-KNOWN-SMLIB-MISSES.md` is mandatory alongside parity. In particular:

- Tiny hill loop and isolated hilltop contact: a matching empty donor result cannot establish emptiness. The independent plane proof reports loop recovery / named tangent refusal.
- Near-parallel plane and box crossings, small/collapsed line intervals and reflected pseudo-box bases: preserve the public negative regressions rather than treating donor pruning as truth.
- Requested-domain endpoint and surface-side answers: retain clipping and typed evaluator failures. Grazing sub-domain regressions were added without a new C++ partial-domain call; that distinction stays explicit.
- Quartic knot-sample iso-line shortcut: require the whole-span control-point bound; successful saved iso rows do not cover this negative case.
- Mixed-curvature saddle and singular start direction: preserve the singular classification and inside branch regression.

The known-miss list also retains invalid construction and absent-solution failures and other typed omissions (curve self-intersection, certain circle rebuilding and undefined zero-tangent paths). The three rows above are not a complete list of all refusals. Historical probe faults are not geometric misses and were not reproduced here. Candidate parity, ULP comparisons and worker-bit agreement cannot replace Catallax acceptance proof.

## Missing decision inputs and completion boundary

Missing: actual B01–B03 specifications; public per-request call attribution connecting the three omissions to named Phase 1 targets; complete pinned-oracle positive outputs for the two refusal rows; and a dedicated seam-start end-to-end case with branch evidence. Public README provenance identifies private raw data but does not make it available or authorized. Firstmate can supply public tracked outputs or commission an authorized local oracle wave; P04 does not access that data or run C++.

This is a bounded documentation/evidence report. No algorithm, Result/build adapter, PRD, tolerance, allowance, fixture expectation or saved/proven answer changed. No qualification of Phase 1.1–1.5, all-thread behavior, 20M work, wall time, import census, edit chain, application workflow or independent STEP import is claimed. The only newly executed behavioral evidence is the Rust saved-gate subset recorded in validation.
