## 2026-10-19: increase the piece 1 default work budget

Pieces: 1
Change: `BudgetPolicy::default().max_work` increases from `100_000_000` to `1_000_000_000`. `max_precision_bits` remains `16384`.
Why: the architect approved the retune after filtered Descartes isolation exhausted the old limit on degree-64 Chebyshev polynomials. Speed is measured, not gated; the work limit still bounds each decision.
Migration: callers using the default receive the larger work budget. Explicit policies are unchanged.

