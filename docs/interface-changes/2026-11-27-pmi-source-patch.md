## 2026-11-27: PMI edits in a STEP file's own records

Pieces: 8; SolidDraft adapter
Change: `catallax_import` exposes `patch_step_pmi` and adds `StepPmi::sources`. The record layer exposes `catallax_step::pmi::patch::{Records, Patch}`. See the [source patch contract](../pieces/08-step-import-repair.md#editing-manufacturing-information-in-place).
Why: ROADMAP step 5.6, plan row C-4D. A narrow PMI edit that keeps unrelated AP242 records, including assembly, colour and opaque metadata.
Migration: Add `sources` to `StepPmi` literals.
