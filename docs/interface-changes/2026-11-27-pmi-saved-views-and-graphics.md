## 2026-11-27: saved views with callouts, and graphical-only PMI

Pieces: 7, 8; SolidDraft adapter
Change: `catallax_step::pmi::Pmi` gains `views: Vec<SavedView>` and `catallax_import::StepPmi` gains `graphics: Vec<GraphicCallout>`. See the [writer](../pieces/07-step-export.md#manufacturing-information) and [reader](../pieces/08-step-import-repair.md#manufacturing-information) contracts.
Why: ROADMAP step 5.6, plan row C-4C. Generalizes the single finish callout to saved-view callouts and reads graphical-only PMI.
Migration: a struct literal of `Pmi` adds `views: Vec::new()`, and one of `StepPmi` adds `graphics: Vec::new()`. S-4I converts SolidDraft's orbit `SavedCamera` (target, distance, yaw, pitch, field of view) into `Camera` and the view's visible PMI IDs into item names.
