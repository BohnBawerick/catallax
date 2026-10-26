## 2026-10-24: native analytic draft sampling

Pieces: features, field, compare

Change: `Plan::visit_feature_candidates` also visits exact rational points on surviving draft facets, cone walls and cap disks. The points come from the field definition, not B-rep vertices or display meshes. Copy roles resolve to the same source patch. Algebraic placements without rational candidates retain `field.sample-algebraic-placement`.

Draft support construction preserves material sense against the source radial normal for oblique planar neutral planes. A neutral orientation without that proof returns the legal unfinished-coverage refusal `draft.unsupported-neutral-orientation`. Exact rational reduction removes common algebraic scales from plane normals, unit directions, cone slopes and apex expressions; neither the surface sets nor budgets change.

Why: comparison requires independent reverse samples on the replacement patches. Unreduced but rational support coefficients unnecessarily exhaust later exact predicates.

Migration: no signature, format, tolerance or work-allowance changes. Candidate visitors may now receive draft points through the existing callback and early-stop contract.

