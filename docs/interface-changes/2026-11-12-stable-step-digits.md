## 2026-11-12: STEP decimals settle before they are written

Pieces: 7
Change: A decimal written by STEP export is accepted only when both ends of its enclosure round to the same text, or an equality proof makes it exact. Otherwise the enclosure narrows before the text is written; a stalled enclosure refuses `step.budget-exhausted`.
Why: Under a smaller precision limit a value next to a rounding boundary could write another last digit, so the same question wrote different STEP bytes at different limits. The same-answer rule needs the written digit to be independent of the limit.
Migration: No caller signature changes. A saved answer may differ in the digits its unsettled enclosure wrote before; this change admits the 13 enumerated replay answers, whose topology and properties are unchanged and whose changed numbers stay inside the answer's stated bound. A value that cannot settle within the limit now refuses instead of writing an arbitrary last digit.
