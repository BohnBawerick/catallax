## 2026-11-27: saved views with callouts, and graphical-only PMI

Pieces: 7, 8; SolidDraft adapter
Change: `catallax_step::pmi::Pmi` gains the field `views: Vec<SavedView>`; a `SavedView` has a name, a `Camera` (eye, direction, up, distance, window, perspective) and the items it shows. Each shown item gets a callout associated with its semantic entity, and each view is a draughting model with a camera. `catallax_import::StepPmi` gains `graphics: Vec<GraphicCallout>` for callouts that present no read item. `read_step_pmi` fills both.
Why: ROADMAP step 5.6, plan row C-4C. Generalizes the single finish callout to saved-view callouts and reads graphical-only PMI.
Migration: a struct literal of `Pmi` adds `views: Vec::new()`, and one of `StepPmi` adds `graphics: Vec::new()`. A file without views writes the bytes it did before. S-4I converts SolidDraft's orbit `SavedCamera` (target, distance, yaw, pitch, field of view) into `Camera` and the view's visible PMI IDs into item names.
