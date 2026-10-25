## 2026-10-24: current-definition boundary candidate handles

Pieces: field, comparison

Change: replace `Plan::visit_feature_candidates` with `Plan::visit_candidates`. Its callback receives `(RPoint3, FaceId)` and covers current primitive, sweep and draft definitions. `Plan::result_face(index)` retrieves a current registry candidate; history and unused entries return `None`. A shared origin may have multiple candidates on different supports. Comparison preserves these handles for reverse samples and tries every current support for a named forward witness.

The obsolete recipe-primitive sampler is removed. Modified faces no longer inherit a historical cutter grid that can place probes on their new sharp edges. Draft trim and sampling predicates match complete origins, not just identically spelled roles from different creators. The native primitive charts retain exact rational positions and share their angular samples within one visit.

Why: a drafted neighbour can split while keeping its origin, and only one piece may move with a merged alias. A single name-to-surface representative loses valid coverage. Same-role coplanar patches also need separate trimmed domains.

Prior art: [OCCT shape history](https://occt3d.com/dev/doc/refman/html/class_b_rep_tools___history.html) distinguishes a source's ancestry from the list of geometric shapes modified from it. Candidate handles preserve that one-to-many distinction without making new recipe names. No source code is copied.

Migration: candidate visitors query the supplied handle directly and read its origin from `plan.faces()[id.index()]`. The old callback and sampler are removed, not retained as aliases. All repository callers are migrated. Primitive frame failures now use the field sampling error namespace. Numeric proof policy, cooperative stop, comparison accuracy and coverage requirements are unchanged.

