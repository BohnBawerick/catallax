## 2026-11-20: semantic PMI in STEP export and import

Pieces: 7, 8; SolidDraft adapter
Change: `catallax_step::pmi` adds `Pmi<F>` and its item types, generic over the face key, with `Pmi::validate` and `pmi::nominals`. `catallax_step::write_step_with_pmi` writes one body with its manufacturing information as AP242 semantic PMI and names every `ADVANCED_FACE`; `catallax_import::read_step_pmi` reads it back into `StepPmi`, keyed by face entity name. New refusal codes: `step.pmi-invalid`, `step.pmi-not-measurable`, `step.pmi-face-not-written`, `step.pmi-unsupported-body`.
Why: the Checkpoint 2 drawings and PMI item ([decisions, section 4.8](../checkpoints/checkpoint-2-decisions.md#48-drawings-and-pmi-the-p1-bracket-test)) needs datums, dimensions and tolerances to travel in the STEP file, attached to faces, so a drawing can be made from the file alone.
Migration: none for existing callers. `write_step` and every other writer keep their bytes. A caller with manufacturing information passes model faces to `write_step_with_pmi` and maps `read_step_pmi`'s face keys through its import's face names.
