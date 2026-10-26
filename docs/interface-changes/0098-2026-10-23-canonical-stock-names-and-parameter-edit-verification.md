## 2026-10-23: canonical stock names and parameter-edit verification

Pieces: recipe, topology, features, field, verification CLI and target scoreboard

Change: `CanonicalStock` supplies exact stock-to-sketch/sweep definitions and cone-cap origin correspondence. `evaluate_canonical_stock` and `Plan::compile_canonical_stock` expose equivalent named results without replacing the fast stock builders. `Model::faces_with_origin` and `edges_with_name` resolve complete result origin sets, including splits and merged contributions. `map_face_origins` rewrites unpublished build provenance without changing geometry or incidence. Field domain predicates and candidate sampling retain generated feature roles independently of renamed public origins.

The result resolvers accept complete nested pattern/mirror `copy` origins. Field feature sampling checks result ownership as well as the defining domain role, so unused seed and mirror-input aliases cannot become result sample names.

`CanonicalStock::new` checks the stock dimensional domain before expansion and returns `Result<Option<CanonicalStock>, Refusal>`. Impossible dimensions report `recipe.impossible-stock`, rather than becoming a legal negative sweep or a reflected meridian. Direct stock and authored feature validation are unchanged.

The opt-in `name-edit` adapter accepts `catallax-name-edit/1` with original and edited recipes. It rebuilds both kernels and reports surviving result-face origins and physical-edge references, missing names, unnamed entities and feature waits. Field registry lookup is labelled separately from boundary proof. `parts/edit-checks.json` specifies P1, P2 and P4 edits; scoreboard rows add `edit_check`, and failed or invalid edits increment `edit_failures`.

Configured checks persist in a top-level `edit_checks` prepass before display or sampled comparison. Report `status` distinguishes running, completed and interrupted coverage. A budget refusal retains its first blocked operation without increasing the allowance. Geometry rows still carry their individual `edit_check`.

Matched-pair timing artifacts record `load_average_before` and `load_average_after` on each pair, beside the pinned-core and alternating-order configuration.

Why: persistent references depend on definitions and incidence, not topology arena positions, dimensions or evaluation-tape handles.

Migration: existing recipes, question adapters, stock geometry, budgets, tolerances and STEP options are unchanged. Generic authored revolve side names are unchanged; cone-cap remapping belongs only to explicit canonical stock expansion. Clients of the scoreboard can inspect `edit_check.status` rather than inferring edit survival from successful construction.

