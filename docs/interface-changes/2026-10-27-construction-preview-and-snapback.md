## 2026-10-27: construction-first preview and snapback

Pieces: application entry, field display and revision lifecycle

The existing `catallax_app::lifecycle::Controller` now supports a construction display before exact B-rep evaluation. It retains the last proved revision independently of the current input. There is still one revision authority and one reserved preview/proof allowance.

### Construction and proof

`Revision::construct(recipe, construction_options, proof_options, allocation, decider)` compiles the immutable definitions with the existing field compiler and requests its unproved display. It does not call B-rep evaluation first. `ConstructionOptions` names visual sampling and resource controls, not certified deviation bounds.

Unsupported definitions return their typed refusal. A preview containing unresolved cells returns `application.preview-incomplete`; it is not published as a complete selectable shape. This route does not substitute an earlier body, fill gaps or silently switch to exact B-rep construction.

The returned `ProofJob` retains the immutable recipe and the reserved proof allowance. When scheduled, it evaluates the exact B-rep through `build_incremental`, reusing the controller's last proved checkpoint only where the existing evaluator establishes an unchanged prefix. Topology and certified display still run before confirmation. No construction is repeated after that retained result is available.

`Revision::preview` remains the route for constructing an exact retained checkpoint before sampling its B-rep preview. Its deferred job meshes that checkpoint rather than rebuilding it. Synchronous `build`, `build_incremental`, `build_soliddraft` and `build_certified` keep their contracts.

### Delivery, cancellation and output gates

- `Controller::latest()` describes only the current delivery. It is absent while construction has not delivered a preview, including an early construction refusal.
- `Controller::displayed()` returns the current unproved/confirmed snapshot, or the last proved snapshot when current proof failed or construction has not delivered. Its revision remains the revision actually displayed. The current failure remains available through `latest()`; an early refusal remains the construction call's error.
- `Controller::confirmed()` returns a body only for the current confirmed delivery. It never supplies the previous successful body for pending or failed input.
- `Controller::cancel()` invalidates the requested edit and restores the last proved snapshot, or an empty state if none exists. The application must also restore the corresponding requested definitions when cancelling its edit.
- A late success or refusal cannot publish after a newer input or cancellation. An unstarted stale job performs no proof work.

Confirmation carries `Arc<IncrementalApplicationBody>`, including its retained checkpoint, names and certified display. Read exact definitions through `body.checkpoint.built()`. This is the existing model/display guarantee, not a universal measurement, whole-boundary, comparison or output certificate. Requested exports and measurements still perform the proof applicable to their scope. An output-only refusal does not change a sound model's status.

Both phase caps are reserved before construction. The unchanged application total is 20,000,000 work units. Unused units do not move between phases or return after cancellation. This entry neither resets an active caller allowance nor adds a second default allowance. `WorkReport` distinguishes reserved caps from arithmetic actually performed.

### Selection migration

`Snapshot::preview()` now returns `Display`:

- `Display::Brep` contains the existing B-rep preview and names. `brep()` borrows those data. Face, edge and vertex selection methods retain their same-revision resolution rules.
- `Display::Construction` contains a `ConstructionDisplay`. `construction()` borrows its mesh and displayed-face origins. Its registry indices are not B-rep handles. Only origins of faces that emitted triangles are offered for selection.

Use `Snapshot::select_origin` to retain a construction-face selection. Resolution waits for confirmation of that same revision and follows its origin against the exact model. A selection cannot cross controllers or revisions, even if numeric identifiers match. Construction previews do not invent exact edge or vertex handles.

### Exercised behavior and limits

The lifecycle regression suite covers same-revision origin resolution, conserved allocations, exact phase boundaries, stale success/refusal, cancellation, incomplete construction refusal, and last-proved snapback without stale export permission.

An API smoke constructed changed block displays before their exact proof, rejected queued old success and failure, cancelled a pending revision, and restored the last proved revision after current proof failure. A guarded caller produced no pending or failed STEP file, then wrote a 6,369-byte STEP after restoring the proved revision. This exercises the gate and actual writer, not independent STEP read-back qualification.

A separate construction-only probe displayed the complete public P2 and P1 definitions with no unresolved cells at an explicit 0.5 mm visual sampling target. It did not prove either exact body. In particular, a full-P1 construction display is not full-P1 B-rep completion. Those shared-host observations do not establish an application click budget.

This kernel entry does not start worker threads, integrate SolidDraft controls or make browser proof responsive. U5 requires the application's revision, selection, cancellation, failure presentation and output paths to adopt the contract and pass their real interaction checks. The earlier exact-build-first measurements remain historical evidence in the [26 October lifecycle record](2026-10-26-application-preview-lifecycle.md).
