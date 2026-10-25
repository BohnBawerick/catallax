## 2026-10-25: evaluated prefixes and indexed construction-fact closure

Pieces: numbers, features and application

Change: `evaluate_incremental` accepts an owned `EvaluationCheckpoint` and
evaluates only the suffix after its unchanged steps. The checkpoint exposes an
immutable `Built` view and retains profiles, body bounds and feature definitions.
Step, format or numeric-policy mismatches rebuild from the start. Application
`build_incremental` adds the existing topology check, display and picking names
without copying the checkpoint's model. Existing build signatures remain valid.
The incremental entry also preserves the existing aggregate application work
guard and any explicit caller-owned allowance.

Why: replaying every earlier boolean made edit cost grow with the feature tree.
Separately, negative fact queries repeatedly traversed equivalent historical
geometry even when the live operands did not change.

Change: `FactBook` indexes equivalence classes and directed incidence between
coincident classes. Union rewires both directions, preserving cycles and the
requirement for at least one actual incidence edge. All recorded facts, geometry
handles and provenance remain unchanged.

Migration: retain a checkpoint per evaluated body and consume it on the next
append. Retain or explicitly clone earlier checkpoints for edits and undo.
Only newly executed suffix work is charged to that request; cached proofs were
paid for when constructed. Full evaluation, certified evaluation, recipe and
manifest formats, STEP options, tolerances and budgets are unchanged. Public
repeated-cut and distinct-hole benchmarks contain timing, strict-budget checks,
per-prefix STEP comparisons and the same-facts query counterfactual.

