## 2026-11-28: assembly PMI without one part is an item diagnostic

Pieces: 8; SolidDraft adapter
Change: `read_step_assembly` replaced the file-level `PMI without a unique part attachment` refusal with item diagnostics. See the [STEP reader contract](../pieces/08-step-import-repair.md#manufacturing-information).
Why: plan row C-2B and the identity rule: refuse malformed controls as items and keep unrelated valid PMI.
Migration: a caller that relied on the refusal reads the parts' `skipped` lists instead.
