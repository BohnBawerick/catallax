## 2026-10-26: application preview lifecycle and budgets

Pieces: application entry, numbers and decisions

Implements the application lifecycle lane in [the preview and deferred proof decision](../ARCHITECTURE.md#preview-and-deferred-proof-decision-of-2026-10-26), using [the unproved preview mesh entry](2026-10-26-unproved-preview-mesh.md). Existing synchronous callers do not migrate. `build`, `build_incremental`, `build_soliddraft`, `build_certified`, their error separation and their 20M default remain unchanged.

### Delivery and scheduling

`catallax_app::lifecycle::Controller::begin()` invalidates the previous pairing immediately and returns a single-use `Revision`. `Revision::preview(recipe, preview_options, proof_options, allocation, decider)` constructs an exact retained build through `build_preview`, publishes a selectable `Status::Unproved` snapshot and returns an owned `ProofJob`. It never calls certified meshing or starts a thread. A revision already invalidated by another edit returns `Ok(None)` without publication. A construction or preview refusal returns `Err` without publishing a partial preview. Callers can keep an old snapshot visible, but must label it with its old revision.

`Controller::latest()` returns an immutable snapshot with its revision, unproved preview mesh and names, status and work report. Geometry and mesh buffers are shared rather than copied on reads. Snapshot fields are private. Numeric handles belong only to the mesh or body in that snapshot.

The caller schedules `ProofJob::prove(decider)` when it chooses. It uses the retained build, not a reread or reevaluation of a mutable recipe. `ProofCompletion::publish()` atomically replaces the status with `Confirmed(Arc<ApplicationBody>)` or `Failed(Error)`. It returns `false` for stale work. A queued old success or refusal cannot replace a new edit, even if that edit has not finished constructing its preview. Jobs and completions cannot be cloned or replayed. Dropping an unstarted job leaves its preview unproved, never confirmed.

`Confirmed` has the existing synchronous `build` guarantees. Construction-time exact feature proofs are unchanged. The deferred phase checks topology and produces a certified display, which the caller can use to replace the unproved mesh. It does not certify the preview triangles, whole-model deviation, a field result or second-kernel agreement. The existing certified inward-offset construction route remains `build_certified`; deferring its construction proofs is a different lane.

A failed display proof keeps its original typed `Error::Display(Refusal)` with the failed preview. No confirmed body is available on failure. Requested measurements and exports use only the confirmed body and retain their own scoped proof and refusal requirements. A measurement- or export-specific refusal does not change this model/display status.

### Selections

`Snapshot::select_face`, `select_edge` and `select_vertex` turn a pick from that snapshot into an opaque, revision-bound `Selection` containing origins or an incident-origin name. Unnamed edges or vertices have no saved selection. `Snapshot::resolve` returns confirmed-result face, edge or vertex handles only after confirmation of that same revision. It resolves names against the proven result, not by copying preview handles. Another edit or another controller rejects the selection even when numeric revision and entity IDs coincide. An empty resolved list means the name has no matching piece. Returned handles remain local to that confirmed body.

### Conserved work allocation

Callers choose explicit construction/preview and proof caps in `WorkAllocation`. `WorkAllocation::application(preview_cap)` splits the unchanged 20,000,000 total; it does not choose an implicit ratio. Public fields let a caller explicitly allocate a different caller-owned total. Caps never borrow unused units from the other phase. Overflow returns `application.invalid-allocation` through the new `Error::Allocation(Refusal)` variant.

Both caps are reserved together before preview construction. An active enclosing operation allowance pays for the full reservation immediately, including the deferred proof cap. Insufficient caller work returns `numbers.operation-budget` through `Error::Allocation`; it cannot fall back to a fresh application budget. Unused or cancelled reserved units are not refunded. This conservative reservation prevents another operation from spending the deferred allowance while proof is queued or after its originating scope returns. Reservation is distinct from actual arithmetic performed. `WorkReport` exposes the allocation and actual `OperationWork` for preview and, after completion, proof. Actual work in either phase cannot exceed its cap. Exhausted phases never publish a confirmed result.

`catallax_numbers::OperationAllowance` supplies this single-use reservation. `reserve(limit)` debits an active operation or owns the explicit standalone limit. `split_off(limit)` transfers units without duplication. Consuming `run` reports actual work under that cap without charging already-reserved units again. It restores any surrounding operation state and caches on normal return or panic. Per-request certificate caches are not transferred between phases. Existing `with_operation_budget` nesting rules and callers remain unchanged.

The caller owns scheduling, decider policies and deadlines. No wall-time measurement affects a geometric verdict. Field preview, certificate reuse, SolidDraft integration and print are unchanged.

### Development measurement

`cargo bench -p catallax-app --bench preview_lifecycle -- "$RESULTS" 3` exercises this lifecycle on the same six fixtures as `preview_mesh`. Set `RESULTS` outside the checkout and pin the command to an available core when comparing samples. Each fixture has one warmup and three measured runs with an explicit 10M preview plus 10M proof split. `preview_wait_ms` includes construction, unproved mesh, edge polylines and selection names through publication. `proof_phase_ms` measures the deferred topology/certified-display job. `proof_completion_ms` measures elapsed time from starting the edit through proof publication, with no artificial scheduling delay. Rows retain phase work, output sizes, typed refusals and exhaustion. These are shared-machine development measurements, not a qualification of the 100/500/1000 ms click or 2000 ms background targets. Raw evidence retained privately.

Observed on 26 October, pinned to one available core, release build, three
measured runs after one warmup. Medians are independent per column. The
construction and proof phases share one 20M request, split 10M/10M.

| Fixture | Selectable preview wait, ms | Deferred proof phase, ms | Proof publication from edit start, ms | Outcome |
| --- | ---: | ---: | ---: | --- |
| fillet | 359 | 2,913 | 3,271 | confirmed |
| chamfer-P2 | 308 | 1,095 | 1,403 | confirmed |
| chamfer | 798 | 1,300 | 2,098 | confirmed |
| countersink | 685 | 1,121 | 1,786 | confirmed |
| P3 | 562 | 3,578 | 4,092 | confirmed |
| P1 | no preview, refusal at 2,723 | not run | not run | `numbers.operation-budget` at `slot-cut-2` |

All five delivered previews preceded proof completion in every measured run.

| Fixture | Actual preview work | Actual proof work |
| --- | ---: | ---: |
| fillet | 1,850,742 | 2,003,783 |
| chamfer-P2 | 3,194,127 | 4,145,970 |
| chamfer | 9,037,463 | 4,152,175 |
| countersink | 2,442,767 | 2,363,169 |
| P3 | 1,852,890 | 1,944,776 |

P1 exhausted the selected preview cap during construction and is not counted as
a delivered preview or a built part. Fillet and P3 proof exceeded the 2,000 ms
background target. These observations do not close their performance gates.
