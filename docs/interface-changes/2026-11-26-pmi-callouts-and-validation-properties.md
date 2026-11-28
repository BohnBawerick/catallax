## 2026-11-26: graphic surface texture callouts and PMI validation properties

Pieces: 7; SolidDraft adapter
Change: `write_step_with_pmi` added surface texture drawings and validation properties. See the [STEP writer contract](../pieces/07-step-export.md#manufacturing-information) for the current encoding.
Why: Checkpoint 3 build order item 1 and research row 27.
Migration: SolidDraft's Open CASCADE expectation for the P1 bracket without saved views gains one `DimensionPresentation` on the bore when it moves its Catallax pin past this change.
