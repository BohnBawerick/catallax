## 2026-10-19: STEP export refuses when the budget allows no bound

Pieces: 7
Change: new refusal code `step.budget-exhausted`. The export refuses with it, writing nothing, when the budget allows no proven bound on a value or spline fit the file needs, or no enclosure narrow enough to write a number within one unit in its last written place. A spline fit with a proven bound that misses `curve_tolerance` is still written with the mark `step.tolerance-not-reached`.
Why: the description gave no outcome for an export whose budget reaches no bound at all, and the export must not invent a `numbers.*` refusal. Piece 6 already has the matching `measure.budget-exhausted`.
Migration: the piece 7 creator returns the code in those cases. The piece 7 test writer accepts it in case 18 and its zero-work variant.

