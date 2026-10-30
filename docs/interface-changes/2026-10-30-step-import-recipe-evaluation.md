# STEP import recipe evaluation

V9-K connects `Operation::ImportStep` to the existing reader and repair implementation. It does not change STEP parsing, repair policy, or STEP writing.

## Sources and builds

`catallax_features::ImportFiles<'a>` maps a source SHA-256 to borrowed immutable bytes. Names in `FileRef` remain labels. The evaluator does not read paths or write source files.

`evaluate_with_files` and `evaluate_incremental_with_files` accept this map. The latter returns the same `EvaluationCheckpoint` as ordinary incremental evaluation. Sources for every import in the recipe must be present and match their saved SHA-256, including imports in a reused prefix. Missing or mismatched bytes return `BuildError::Recipe` with code `recipe.parse` and a step-local explanation.

`catallax_app::build_incremental_with_files` accepts the map and returns an `IncrementalApplicationBody`, with exact geometry, display, and picking names. It uses the same bounds, work allowance, and checkpoint rules as `build_incremental`. Existing file-free entries retain their signatures. An import through those entries fails as a missing-source recipe error.

The evaluator calls `catallax_import::apply` with default import options and the recipe step ID. It places the repaired body using the step's pose. Face origins retain `{step, imported: {key}}`, where each key is the source `ADVANCED_FACE` entity name. A measured enclosing box supplies the conservative radius used by later feature operations.

## Saved decisions and refusals

Only `decline` is supported in the saved decisions map for this batch. Keys use the kernel form `repair:<first 12 source hash digits>:g<index>`. An automatic repair saves no answer. Explicit acceptance, unknown questions, and questions about other bytes refuse with `import.invalid-decision`.

A step must yield exactly one solid. `import.multibody` reports the count and source solid entity names, including any left open by a decline. `import.left-open` names the declined solid and returns no body. Reader and repair failures retain distinct `import.*` codes and their entity details. No failed evaluation publishes a partial checkpoint.

Repair records remain the existing `catallax_import::apply` contract. This build entry does not add them to the recipe or change their persistence format.

## Acceptance and companion

`crates/catallax-app/tests/import_step.rs` covers the independently specified K1-K3 cases on the committed V9 fixtures. Replay R remains the full 477 earlier answers in one worker configuration, with direct no-healing STEP checks and the approved known-raw-invalid list.

SolidDraft V9-S must pin the merged kernel, supply source bytes, and complete preview, saved decisions, reopen, and checked export. Its question-id helper currently emits `repair:<hash>:0`, while the kernel emits `repair:<hash>:g0`. That mismatch is in SolidDraft and is not changed here. V9-K alone does not complete the V9 user scenario.
