## 2026-10-22: named B-rep sketch features and recipe version 3

Pieces: recipe, topology, boolean, features, measurement; the judge exchange is unchanged

Change: add `catallax-features` with exact sketch validation, odd-nested profiles, signed extrusion, analytic revolution, recipe stock-solid inputs and explicit booleans. `Operation::ExtrudeThrough` defines a two-sided through-all cut whose finite construction caps are proved outside the target. The schema becomes `catallax-recipe/3`. `Model::add_named_face`, `face_origins`, `edge_name` and `vertex_name` expose canonical recipe provenance without changing `Face` layout or unnamed legacy builders. Boolean splits retain origin sets; maximal face merges union them.

Builds retain an immutable recipe snapshot so operation definitions remain available after analytic reduction. Full revolutions reject isolated axis contacts that would pinch a shell; valid partial sectors remain buildable. Straight cone generators ending at an apex contribute exact zero to the boundary integrals, avoiding a removable singularity without raising measurement budgets.

`catallax_boolean::make_maximal` reuses the sewing merger on constructed boundaries while preserving their proven shell/lump nesting and singular loops. Half-turn caps joined along the axis merge with both origins; artificial axis edges and degree-two vertices disappear. The existing boolean stages retain their previous classification and refusal behaviour.

Why: recipes need actual B-rep features and persistent operation names. Rational analytic representations are an explicit intermediate implementation boundary, not new impossible-input rules. RECIPE.md section 7.14 lists the remaining legal-input bugs and their specific unsupported-geometry codes.

Migration: explicitly relabel reviewed version 1 or 2 recipes as version 3, then parse/check/write. Earlier operation meanings are unchanged. Match the new operation variant; field evaluation of sketch features and through-all remains unsupported. Use `catallax_features::evaluate` for named B-rep recipe builds, then the existing measurement and STEP APIs on its returned model and body. No command, question format, manifest format, legacy unnamed geometry or STEP format changes.

