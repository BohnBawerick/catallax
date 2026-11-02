# Planar import local tolerances

Pieces: STEP import, feature evaluation, application checkpoint builds.

## Contract

`catallax-import-repair/2` replaces `/1`. `RepairRecord::face_tolerances` is required, including when empty. Each `FaceTolerance` names the immutable source solid and face, encloses the whole trimmed-face two-sided Hausdorff deviation with `deviation_lower_mm` and `deviation_upper_mm`, and records `approved_tolerance_mm` and `feature_size_squared_mm2`. Distance endpoints and limits are exact decimal strings in millimetres. Feature size squared is an exact rational string in square millimetres. The attained upper endpoint is the local tolerance, not the requested limit or maximum vertex movement.

The proof covers straight, strictly convex planar polygons with one boundary. It checks their complete covers and uses exact distances to the trimmed polygon in both directions. Nonconvex faces, holes and unsupported nearest-feature covers refuse rather than receive unproved tolerances. Clean imported curved faces retain their existing path. Every changed face must meet its approved file cap and a tenth of its conservative nearest-feature size. User acceptance can increase the file cap, never the feature cap. Final bounds refer to the original source after all repairs; no intermediate tolerance is added or clamped.

The feature bound includes edge lengths and distances of nonincident vertices to supporting edge lines. This catches narrow skew polygons that have long edges. Final certification includes source and repaired neighbouring faces. The existing 0.01 mm automatic threshold and the tenth-of-feature comparison are inclusive.

`Built::import_record(&StepId)` returns the checked source-bound record for an immutable import body. It survives cold builds and checkpoint reuse and remains associated with the original step, not any later modified result. Rigid placement preserves distances. This is not a whole-body material certificate or a certificate for display or output rounding. The STEP writer's existing accuracy contract remains relative to the repaired model and is unchanged.

## Saved decisions

Under each exact `Proposal.question`, recipe decision values now accept:

- `decline`, including declines of automatic proposals.
- `accept:<tolerance-mm>:<approver>` for an asked seam proposal. The tolerance must have an exact decimal representation and cover the measured gap. The approver must be nonempty. Colons after the first approver character remain part of the identity.

An automatic proposal still rejects explicit acceptance. Openings still cannot be accepted. Other-source, unknown and duplicate lower-level decisions retain their refusals. Recipe serialization remains unchanged because decision values were already strings. Replacing source bytes still invalidates the hash-keyed decisions. No app-supplied acceptance is synthesized by the evaluator.

## Migration and SolidDraft

Persist unchanged source bytes, hash, recipe decisions and the `/2` record. Call `rebuild` on reopen to verify every recorded bound and derivative hash. Existing `/1` records are not distance evidence and are refused by `rebuild`; reapply their source and saved approvals or declines and persist the newly checked `/2` record. Do not change only the format tag or invent missing bounds.

SolidDraft must pin this revision, use the exact proposal question, preserve the approver and tolerance in accepted recipe values, and retain `Built::import_record` with the import checkpoint. Show approval limits separately from attained face tolerances and vertex movement. Save declines and report that the affected solid stays open. Exercise its real rendered viewport, save/reopen and checked-export flow before claiming application rollout. No SolidDraft files or its clean analytic rollout are changed by this slice.

## Next slices, in order

1. Version the seam exchange contract and qualify the next periodic analytic imports, including cone-angle and sphere/torus pole handling.
2. Extend curved repairs and nonconvex, holed or spline face covers with complete proofs and source-recorded licensed fixtures.
