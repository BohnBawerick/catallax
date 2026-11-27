## 2026-11-27: face notes and knurl designations as semantic text

Pieces: 7, 8; SolidDraft adapter
Change: `catallax_step::pmi::ItemKind` gains `FaceNote { text, feature }` and `Knurl { designation, feature }`, written as CAx-IF semantic text on their faces and read back by `read_step_pmi`. `ItemKind::is_text` tells them apart from semantic PMI. A semantic text property with an unknown variable name is skipped with its reason.
Why: ROADMAP step 5.6, plan row C-2N. A knurl and a note on the same face must stay separate and keep their faces.
Migration: an exhaustive `match` on `ItemKind` adds the two arms. A file with neither keeps its bytes and its counts. SolidDraft maps its knurl and face note kinds to these variants once S-2N defines them.
