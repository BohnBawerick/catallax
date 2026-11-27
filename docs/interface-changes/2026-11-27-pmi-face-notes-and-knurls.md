## 2026-11-27: face notes and knurl designations as semantic text

Pieces: 7, 8; SolidDraft adapter
Change: `catallax_step::pmi::ItemKind` gains `FaceNote { text, feature }` and `Knurl { designation, feature }`. See the [STEP text encoding](../pieces/07-step-export.md#manufacturing-information).
Why: ROADMAP step 5.6, plan row C-2N. A knurl and a note on the same face must stay separate and keep their faces.
Migration: Add the two arms to exhaustive matches. Use `ItemKind::is_text` to keep the two kinds out of semantic PMI counts.
