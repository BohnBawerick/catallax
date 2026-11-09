# 2026-11-10: rational path schema and loft migration

Pieces: recipe, features, application

This is the schema prerequisite for PLAN stage 5, not loft or spline-sweep construction. Checked rational-patch carriers, their consumers and the operations land separately.

## Saved definitions

Versions 4 and 5 accept `PathSegment::Spline` with fields in this order: `kind`, `id`, `degree`, `control_points`, `weights`, `knots`, optional `corner_axis`. Coordinates, weights and expanded knots remain exact rational strings. Degree is a positive JSON integer. The first control equals the preceding endpoint exactly; the last becomes the endpoint for the following segment. Interior knot multiplicities may equal the degree, but may not exceed it. No internal span gets a new authored id.

The existing loft fields remain `sections`, each containing `profile` and an ordered per-loop `start` list, followed by `shape: "ruled" | "smooth"`. The lists define loop correspondence, independent of stored loop order. Every section must contain the same number of loops, with exactly one start from each loop. Unequal authored segment counts remain valid. No nearest-point matching is inferred.

`Recipe::loft_parameterization()` reports `ArcLength` for version 3 and `RationalCharts` for versions 4 and 5. The latter means the rational charts, exact common refinement and uniform station laws in RECIPE section 4.12. This query does not build or certify a patch.

`Recipe::upgrade_to(format)` validates both versions and returns a new recipe. It never mutates the saved source. Upgrading a version 3 loft, thicken, sweep, fillet, chamfer or draft returns `recipe.migration-review`, conservatively requiring explicit review of the changed definitions. This includes all-line lofts rather than guessing when parametrisations agree. Retain version 3 when the old definition is intended. After review, an editor may explicitly author a new version 4/5 recipe; that is not proof that it defines the old solid. Version 4 to 5 preserves references and operation parameters. Downgrades require review too.

Saved application feature instances transform all spline controls as points. Weights and knots stay unchanged. A corner axis transforms as an axial vector, including the sign reversal under reflection. Sweep orientation remains explicit: `fixed` means pure translation, while `rotation_minimising` means the RECIPE frame law. Neither mode is silently substituted for the other.

## Current construction limits

The feature builder returns `features.unsupported-path` on spline paths, including through the ordinary application entry. It does not approximate the path with a polyline. Loft construction retains `features.unsupported-operation`. Existing line/cardinal-arc sweep reductions are unchanged. The field does not gain spline sweep support. No approximation tolerance, mesh certificate, measurement or export capability is claimed by this schema change.

## Case inventory

Tests are in `crates/catallax-recipe/tests/spline_paths.rs`, `loft_migration.rs`, and `crates/catallax-app/tests/spline_path_schema.rs`.

| Case | Test |
| --- | --- |
| Weighted rational path, repeated interior knot, both saved orientations, versions 4/5 | `weighted_repeated_knot_paths_keep_their_exact_definition_in_both_modes` |
| Version 3 spline rejection | `version_three_rejects_spline_paths` |
| Degree zero, degree at control count, integer overflow, negative degree | `zero_degree`, `degree_at_control_count`, `degree_overflow`, `negative_degree` |
| No controls or wrong knot/weight count | `empty_controls`, `wrong_knot_count`, `wrong_weight_count` |
| Zero/negative weights or JSON numeric weight | `zero_weight`, `negative_weight`, `numeric_weight` |
| Decreasing or out-of-domain knots | `decreasing_knots`, `outside_domain` |
| Unclamped endpoints or excess endpoint multiplicity | `unclamped_start`, `unclamped_end`, `extra_start_knot`, `extra_end_knot` |
| Excess interior multiplicity | `interior_knot_multiplicity_cannot_exceed_the_degree` |
| Exact joins across line/spline/arc segments | `joins_follow_the_exact_last_control_across_mixed_segment_kinds` |
| First control differs by a tiny exact amount | `first_control_must_equal_path_start_without_a_join_tolerance` |
| Repeated controls and stationary points | `repeated_controls_and_stationary_points_are_not_schema_errors` |
| Optional corner axis and duplicate authored id | `spline_corner_axes_are_retained_and_ids_stay_unique` |
| Unknown spline field | `unknown_spline_field` |
| Unequal segment counts, authored loop order, ruled/smooth laws | `rational_lofts_accept_unequal_segment_counts_and_keep_authored_loop_order` |
| Unequal loop counts | `loop_counts_must_agree_between_sections` |
| Missing, repeated or wrong-section start | `a_missing_loop_start_is_not_a_loop_count_error`, `repeated_start_in_one_loop_does_not_cover_another_loop`, `a_start_from_another_section_is_not_a_correspondence_guess` |
| Fewer than two loft sections | `a_loft_requires_at_least_two_sections` |
| Legacy loft law survives read/write | `legacy_loft_reads_and_writes_without_changing_its_law` |
| Legacy loft upgrade is not automatic | `legacy_loft_cannot_be_automatically_upgraded_even_with_only_lines` |
| Version-dependent feature defaults and sweep modes require review | `version_dependent_feature_defaults_also_require_review`, `legacy_sweep_transport_requires_review_in_both_orientation_modes` |
| Unchanged graph upgrade and version 4 to 5 loft upgrade | `unchanged_operation_graph_can_upgrade_without_rewriting_names`, `reviewed_rational_loft_upgrade_keeps_stations_starts_and_shape` |
| Downgrade, unknown format and invalid source | `migration_rejects_downgrades_unknown_formats_and_invalid_input` |
| Reflected and translated saved spline definition | `reflected_translated_path_controls_and_half_turn_axis_save_exactly` |
| Application refuses rather than sampling an unsupported path | `spline_paths_are_not_replaced_with_sampled_line_sweeps` |

## SolidDraft follow-on

No SolidDraft files change here. Its future Loft command must save ordered section references, explicit per-loop start ids and ruled/smooth selection in version 4 or 5. Its Sweep command must save spline controls, positive weights, expanded knots and the chosen orientation. Legacy lofts must keep version 3 unless their correspondence is explicitly reviewed. Kernel schema acceptance must not enable either command before the checked construction, display, naming and requested-measurement paths land.
