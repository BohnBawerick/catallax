## 2026-11-28: assembly PMI without one part is an item diagnostic

Pieces: 8; SolidDraft adapter
Change: `read_step_assembly` no longer refuses the whole file with `Unsupported` and the type name `PMI without a unique part attachment`. An item whose links reach no part, or more than one, is listed in each part's `StepPmi::skipped` with its reason, and an item whose links reach one part is read, or refused with its reason, there. A shape definition of a property, such as a broken datum target's parameters, is never taken for a product's geometry.
Why: plan row C-2B and the identity rule: refuse malformed controls as items and keep unrelated valid PMI.
Migration: a caller that relied on the refusal reads the parts' `skipped` lists instead.
