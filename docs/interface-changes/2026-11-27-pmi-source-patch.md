## 2026-11-27: PMI edits in a STEP file's own records

Pieces: 8; SolidDraft adapter
Change: `catallax_import::patch_step_pmi(bytes, &[PmiEdit])` replaces a geometric tolerance value in the source records and keeps every other record byte for byte, refusing with the code `step.pmi-patch` before any output when the edit is not proven safe. `StepPmi` gains `sources`, the entity each item was read from. `catallax_step::pmi::patch` exposes the record layer: `Records` and `Patch`.
Why: ROADMAP step 5.6, plan row C-4D. A narrow PMI edit that keeps unrelated AP242 records, including assembly, colour and opaque metadata.
Migration: a struct literal of `StepPmi` adds `sources`.
