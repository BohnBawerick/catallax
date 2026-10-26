## 2026-10-24: benchmark refusal medians and checker provenance

Pieces: verification tooling and milestone data

Change: head-to-head summaries exclude refused phases from successful medians, retain separate phase refusal records and distinguish STEP opening from topology/face-count acceptance. The pinned OpenCASCADE 8 checker uses `TopoDS.Edge`. The first fixed-allowance report retains original raw rows plus `catallax-head-to-head-property-recheck/1` numeric rechecks and SHA-pinned `catallax-head-to-head-summary/1` derived rows.

Why: a fast refusal is not a displayed part, and a checker API error is not a successful correctness check. Correcting independent property checks must not rewrite measured timing or invent STEP evidence.

Migration: regenerate summaries from raw `pairs`, not an older embedded `summary`. Preserve original checker errors beside any later recheck. Kernel interfaces, modelling tolerances and work allowances are unchanged.

