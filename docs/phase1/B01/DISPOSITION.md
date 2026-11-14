# B01 disposition: no demonstrated Phase 1 projected-drop consumer

Close B01 without a port or oracle regeneration for the current Phase 1 scope. Firstmate's B01 assignment explicitly permits this disposition when the merged P04 consumer map and current callers establish no Phase 1 consumer. This is a bounded scope decision, not completion of the projected side-search algorithm.

## Frozen baseline and provenance

Repository: `BohnBawerick/catallax`. Remote main was observed with `git ls-remote origin refs/heads/main`, fetched, and frozen on 2026-11-14 UTC to `e6abb63dcdd1ed49f66ce82065365834684e0c51`. That commit is the merged P04 PR #601, titled “P04: remaining SMLib port calls mapped to Phase 1 consumers.” The fresh isolated branch is `phase1/B01-projected-drop-disposition`. P04's branch, evidence and historical baseline are preserved.

Read the merged [P04 verdict](../P04/PORT-CONSUMERS.md), [known-miss status](../../M9-KNOWN-SMLIB-MISSES.md), [Phase 1 roadmap](../../../ROADMAP.md) and current public Rust callers. No repository `AGENTS.md` or applicable `.agents/skills/SKILL.md` was found in the checkout or shared guidance directory. The existing verification skill was reviewed for context; its CLI qualification workflow is not invoked for this static disposition.

## Current call evidence

Line numbers below refer to the frozen baseline. A repository-wide Rust symbol search for `drop_and_trim_curve` and `coincident_ranges`, followed by inspection of call arguments and wrapper definitions, found:

| Boundary or caller | Inspected evidence | Consequence for B01 |
| --- | --- | --- |
| [Projected-drop boundary](../../../crates/catallax-trace/src/surface_drop_curve.rs), lines 1262–1273 | `drop_and_trim_curve` immediately returns `SolverRefusal::NotPorted(Intersect)` when `on_surface` is false. | The omission is explicit and remains unchanged. |
| [Surface intersector range drop](../../../crates/catallax-trace/src/surface_trace.rs), line 2459 | Drops the side-range curve onto each surface with final arguments `true, true` (`keep_all`, `on_surface`). | Existing side-range consumers do not request the projected false branch. |
| [Curve/surface coincident-range classification](../../../crates/catallax-trace/src/curve_surface_global.rs), line 1226 | Direct call ends in `true, true`. | Does not request projected drop. |
| [Coincident-range wrapper](../../../crates/catallax-trace/src/surface_drop_curve.rs), lines 1320–1329 | `coincident_ranges` always calls the drop helper with `on_surface=true`. Its final own argument is `drop_failure_is_empty`, not `on_surface`. | Calls through this wrapper cannot select the projected false branch. |
| Wrapper callers: [curve/surface](../../../crates/catallax-trace/src/curve_surface_global.rs), line 1340, and [line/surface](../../../crates/catallax-trace/src/line_surface.rs), line 736 | They use `coincident_ranges` after coincidence checks. The line caller's `false, false` arguments are `keep_all` and `drop_failure_is_empty`. | These booleans are not evidence of a projected-drop consumer. |
| [Public drop gate](../../../crates/catallax-trace/tests/drop_curve.rs), lines 68–90; [generator description](../../../crates/catallax-trace/tests/fixtures/drop-curve/generate.py), line 149 | Test inputs can select `on_surface=0`; the bump iso-curve is a `drefuse` row. The harness expects typed refusal for that row. | A test of the omission establishes neither a production consumer nor a successful projected answer. No expectations change. |

The already-ported [general curve drop](../../../crates/catallax-trace/src/surface_drop_general.rs) is a different operation. Its existence, cylinder seam walks and successful saved cases do not establish the missing projected side search. Closed-seam starts and tracer fallback can still refuse inside the true branch; their disposition is outside B01.

## Phase 1 relevance and reason to stop here

Phase 1.3's spline Boolean route uses `Candidates::generate` in [traced pairs](../../../crates/catallax-intersect/src/pair/traced.rs), then `Intersector::do_intersection` in [trace candidates](../../../crates/catallax-intersect/src/pair/trace_candidates.rs). The inspected intersector side-range calls use the true branch above. Phase 1.5 depends on that intersection work, but the public callers supply no projected false-branch demand. No current Phase 1.1 numerical-result or Phase 1.2 import consumer of that false branch was established either. This confirms P04's bounded source verdict; it does not claim every Phase 1 request was executed or all its blockers were identified.

A port now would implement an unneeded branch for the evidenced current consumers. Oracle regeneration is also unnecessary for this no-consumer closure: there is no proposed changed algorithm or positive expectation to compare. The existing refusal stays visible. The [saved drop fixture provenance](../../../crates/catallax-trace/tests/fixtures/drop-curve/README.md) identifies the pinned donor as SMLib `c9c5979f04cea8c38213004d4f04df0dcfb40d05` and says the omitted positive raw output is retained privately. B01 neither accesses that raw evidence nor converts refusal into parity success.

## Concrete trigger to reopen

Reopen B01 when a named authorized Phase 1 recipe or operation supplies a current production call path that needs `drop_and_trim_curve(..., on_surface=false)`, together with a bounded public input, requested curve interval and surface domain, and the required output contract. Confirm why the existing true/general-drop routes cannot meet that consumer. Scope only the needed branch. Obtain missing complete positive donor outputs through the external M/B oracle lane, not cloud regeneration, and retain requested-domain and known-miss controls. Candidate parity must remain separate from geometric acceptance; this disposition authorizes no change to tolerance, allowance, caps or Result/build contracts.

## Validation and limits

Executed only static inspection and document checks for this B01 branch: repository-wide Rust symbol search, call-argument reads, merged P04/baseline history inspection, relative-link existence checks, machine-path audit and `git diff --check`. All linked files exist at the frozen baseline. No tests, builds, toolchain probes, C++ calls, oracle regeneration, dynamic call tracing, Boolean replay, STEP export or qualification were executed. P04's historical Rust test results are not new B01 results.

Only this disposition is changed. PRD, N01 adapters, algorithms, shared design, saved expected answers, other jobs and P02's hold remain untouched. No disposable B01 build output or scratch was created; checkout, caches and required evidence are retained. Publication is a branch push for Firstmate to review; no PR is created here. Exact pushed head and remote verification are supplied in the handoff.
