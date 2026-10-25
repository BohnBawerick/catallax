## 2026-10-22: comparison consumes field profile sweeps

Pieces: field queries, display comparison
Change: add `Plan::visit_feature_candidates` as a query-only, cooperatively stoppable source of rational supporting-face candidates for compiled sketch sweeps. Comparison now accepts sketch, extrusion, through-all and revolution steps in the common rational support subset. Candidate generation propagates conservative field-result bounds through the expression tree to sample through-all walls without finite construction caps.
Why: the field evaluator gained these operations while the comparison library was being built. They must use the same two-way checks as stock solids rather than retaining obsolete unsupported-operation reports.
Migration: call `face_query` before accepting a candidate as a surviving boundary point. The visitor returns `Ok(false)` when its callback stops, not a complete-coverage claim. `field.sample-options`, `field.sample-unbounded`, `field.sample-role` and `field.sample-algebraic-radius` are sampling failures, not geometry refusals or agreement. Existing builders, query meanings, recipe schema and exchange formats are unchanged.

