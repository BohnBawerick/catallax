## 2026-10-18: ND10 properties enclose the built body, never the STEP copy

Pieces: 6, judge interface
Change: `properties` intervals enclose the exact properties of the built body. They contain the true answer's values unless a mark from outside pieces 6 and 7 is present; marks with codes `measure.*` or `step.*` never change the body. This accepts the judge's provisional behaviour for answers without such marks.
Why: a mark from a not-exact decision means the built body may differ from the true answer, so its intervals cannot promise the true values.
Migration: none for piece 6; its enclosures were already of the built body.

