## 2026-10-27: STEP import `rebuild` from source and record

Pieces: `catallax-import`, application adapter (SolidDraft)
Change: Adds `rebuild(model, source, &ImportOptions, &RepairRecord, step) -> Result<Imported, ImportRefusal>`, the reopen path of [`apply`](2026-10-27-step-import-propose-apply.md). The record's approvals become decisions (`automatic` none, `user` accept at the recorded tolerance and approver, `declined` decline), the same deterministic repair runs, and the result must equal the record, derivative hash included. Source bytes whose SHA-256 differs from the record's, or a derivative hash that differs, give `hash_mismatch`; a record whose changes differ from the repair gives `validation_failed`; an unknown `format` gives `malformed`.
Why: an approved change record must survive saving and reopening, with the repaired derivative reproducible from the immutable source and the record alone.
Migration: None. An application persists the record JSON and calls `rebuild` on reopen instead of asking again.
