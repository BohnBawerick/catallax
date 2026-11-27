## 2026-11-27: datum targets and composite or projected position in STEP PMI

Pieces: 7, 8; SolidDraft adapter
Change: `catallax_step::pmi::ItemKind` gains `DatumTarget(DatumTarget<F>)` (point, line, circular or rectangular area on one planar face) and `Position(PositionControl<F>)` (a single frame with a projected zone or datum material boundary, or a composite frame of two segments). `ItemError` gains `DuplicateTarget`. `ItemKind::faces` also lists target and zone entry faces. `write_step_with_pmi` first checks these items against the model with `pmi::check_on_model`, and `read_step_pmi` reads them back.
Why: ROADMAP step 5.6, plan row C-2B. A typed round trip of composite segments, projected zones and face-linked targets.
Migration: an exhaustive `match` on `ItemKind` adds the two arms. A single position with no projected zone and no datum boundary stays `Tolerance` with `Characteristic::Position`; the validator refuses that form as `Position`. A reader of `StepPmi::skipped` no longer sees placed datum targets or datum material boundaries there.
