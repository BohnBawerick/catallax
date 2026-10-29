## 2026-10-27: STEP import and repair wire types and call sequence

Pieces: new crate `catallax-import`, application adapter, recipe `import_step`
Change: Adds the crate `catallax-import` with the serde wire types that cross the application boundary for STEP import, per [piece 8](../pieces/08-step-import-repair.md). This change adds the types only. The functions below arrive in the slices named beside them; until each is merged a caller can code against the JSON but cannot call it. Numbers are decimal strings in millimetres.

Call sequence an application (SolidDraft) follows:

1. `propose(source: &[u8], options: &ImportOptions) -> Result<ProposalSet, ImportRefusal>` (slice 3). Reads the Part 21 bytes, builds the shells, measures gaps and applies the policy table. Returns `catallax-import-proposals/1`: `source_sha256` and `proposals`, each with `id` (`g0`, `g1`, ...), `question`, `entities`, `kind` (`seam_close`, `body_gap`, `opening`), certified `distance_min_mm` and `distance_max_mm`, `tolerance_mm`, `disposition` (`auto` or `ask`) and `reason`. No geometry moves in this call.
2. The application shows every `ask` proposal and collects one `Decision` per proposal: `{"decision":"accept","id","tolerance_mm","approver"}` or `{"decision":"decline","id"}`. `auto` proposals take no decision. An accepted tolerance must be at least `distance_max_mm`.
3. `apply(source: &[u8], name: &str, decisions: &[Decision]) -> Result<Imported, ImportRefusal>` (slice 3). Applies automatic repairs and accepted decisions, validates the result (closed, non-self-intersecting, oriented, connected, measured movement within tolerance) and returns the body, the `RepairRecord` and its derivative hash. A declined proposal leaves the original geometry unchanged and is recorded with approval `declined`.
4. The application saves the source bytes unchanged, the record JSON beside the recipe, and in the recipe step `decisions[question] = "accept:<tolerance_mm>"` or `"decline"`. It never writes an export over the source file.
5. On reopen, `rebuild(source, &RepairRecord)` (slice 4) reproduces the derivative and its hash from the source and record alone. A hash mismatch is `hash_mismatch`, a recipe error.

Record `catallax-import-repair/1`: `format`, `source {name, sha256}`, optional `derivative {sha256}` (absent when every change was declined), and `changes`, each with `id`, `entities`, `kind`, `bounds_before`, `bounds_after`, `max_movement_mm` (measured on the result), `tolerance_mm`, `reason`, `approval` (`automatic`, `user`, `declined`) and optional `approver`.

Typed refusals, JSON `{"code": ...}`: `malformed`, `unsupported {entity, type_name}`, `unknown_unit`, `hash_mismatch`, `not_solid` (geometry), `invalid_decision`, `validation_failed`, `resource_limit` (incomplete, not impossible). Only `not_solid` is a geometry claim.

Why: SolidDraft builds only against merged, documented interfaces, and the import consumer needs one agreed shape before the reader exists.
Migration: None for existing callers. SolidDraft adds the import command and dialog against these types.
