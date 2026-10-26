## 2026-10-23: comparison queries propagate their original refusal

Pieces: compare, verification CLI
Change: `Report` gains `refusal: Option<catallax_numbers::Refusal>`. A failed field face or coactive-face query stops further reverse queries and publishes its original refusal instead of panicking. The verification adapter returns that code and detail; a whole-request numeric-budget refusal still takes precedence.
Why: the first P5 countersink exhausted its unchanged request allowance inside reverse face queries and previously became `verify.internal-error`.
Migration: report constructors initialise `refusal`; adapters propagate a present refusal without converting it to a proven disagreement or internal error. No allowance or tolerance changes.

