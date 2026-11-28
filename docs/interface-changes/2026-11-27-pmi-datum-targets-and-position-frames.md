## 2026-11-27: datum targets and composite or projected position in STEP PMI

Pieces: 7, 8; SolidDraft adapter
Change: `catallax_step::pmi::ItemKind` gains `DatumTarget` and `Position`, `ItemError` gains `DuplicateTarget`, and `PartNotes` gains `standard`. See the [STEP writer contract](../pieces/07-step-export.md#manufacturing-information).
Why: ROADMAP step 5.6, plan row C-2B. A typed round trip of composite segments, projected zones and face-linked targets.
Migration: Add the new arms to exhaustive matches. Add `standard: None` to `PartNotes` literals that use the default standard. Use `ItemKind::faces` when enumerating all attachments, including datum target and projected-zone entry faces.
