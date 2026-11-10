# Rational loft cases

The rational loft entry is `catallax_app::build` with a version 4 or 5 `loft` operation. Its immutable definition is the recipe, including ordered start ids and `shape`. The correspondence law is [RECIPE 4.12](RECIPE.md#412-loft); this inventory follows [PLAN 3.8](PLAN.md).

Tests live in `crates/catallax-app/tests/rational_lofts.rs`. They construct actual recipes through the application entry, not manually assembled `Built` values.

The public successful recipe is [`loft-spline-prism.json`](../crates/catallax-features/tests/fixtures/loft-spline-prism.json). The saved [`loft-natural-circles.json`](../crates/catallax-features/tests/fixtures/loft-natural-circles.json) retains the documented default-budget limit; it is not advertised as a successful eager application build.

| Case | Contract test | Required result |
|---|---|---|
| Two-section ruled and smooth polygon lofts | `two_section_ruled_and_smooth_frustums_have_exact_mass_and_physical_topology` | Both laws give the same frustum; volume 28, centroid height 51/28, six physical faces, twelve edges and eight vertices. |
| Three circular sections, natural interpolation | `natural_three_circle_loft_has_a_documented_application_budget_limit` | Honest `ResourceLimit` at the unchanged 20M application allowance; profiling and completion are a separate follow-up, not a relaxed work cap. |
| Unequal physical station distances | `unequally_spaced_natural_circles_have_a_documented_application_budget_limit` | The same honest default-budget refusal; the follow-up must retain uniform integer station spacing, not substitute chord-length interpolation. |
| Quadratic and cubic spline profiles | `quadratic_and_cubic_profiles_keep_exact_volume_and_area_centroid` | Exact volumes 16 and 144/5; centroid ordinates 4/5 and 9/7; numerical chart cuts introduce no physical entities. |
| Later sections have more authored segments | `later_segment_splits_add_no_names_edges_or_vertices` | Unequal segment counts are accepted; only the first section's four side names survive. |
| First section has collinear authored subdivisions | `first_section_collinear_ids_survive_on_maximal_logical_faces` | Each maximal side retains both contributing first-section ids; no normalization-based reallocation of u. |
| Nonaffine parallel sections | `nonaffine_parallel_sections_use_a_whole_polar_embedding_proof` | A non-parallelogram target is admitted by a whole-domain angular-Jacobian proof; volume 23/3. |
| Nonparallel affine sections | `nonparallel_affine_sections_keep_the_spatial_volume_map` | The complete spatial volume map is embedded; volume 1600/101 for the specified tilted square. |
| Holes and nested islands | `nested_holes_and_islands_remain_separate_material_lumps` | Odd nesting survives; two material lumps, seven faces, six ring edges, no seam vertices, volume 39π. |
| Coplanar sections | `coplanar_sections_regularize_to_empty_without_sheet_handles` | Regularized empty body, zero volume and no selectable sheet entities. |
| Reversed correspondence | `reversed_correspondence_is_a_typed_refusal_not_a_twisted_approximation` | Typed `loft.admissibility-unproven`, with no approximate or partial body. |
| Legacy version 3 arc-length definition | `legacy_arc_length_lofts_are_not_reinterpreted_as_rational_charts` | Typed `loft.legacy-parameterization`; a saved definition is never silently changed. |

Equal-angle circular-arc construction has separate exact unit cases in `catallax-features/src/loft/charts/tests.rs`: one through four subarcs in both directions and cardinal endpoint boundaries. These check exact endpoint rotation and unit-circle identities, including the cubic root needed for three equal subarcs.

## Admission and remaining coverage

Every successful loft has an exact positive-weight rational definition. Its side charts have complete local embedding certificates. G1-adjacent charts form logical faces; exact common-denominator refinement and numerical chart cuts are not selectable edges. The cell topology is retained separately for certified display, integration and STEP face-region mapping.

Whole-body admission is separate from local chart regularity. Parallel affine section families use strict station progression and positive planar determinants. Parallel nonaffine families require a strict whole-domain polar Jacobian and complete inter-loop separation. Nonparallel affine families require a globally strongly monotone projected volume map over a convex box containing the initial material. No sampled intersection test admits a body.

Exact coefficient identities reduce straight translated section families to existing prisms when their profile carriers permit it. Positive rational control-hull proofs reduce planar logical sides and straight physical boundaries to ordinary planes and lines. These are exact reductions of the recipe definition, not fits, and avoid running rational-surface integration and display on analytic geometry.

These certificates are sufficient conditions, not a general loft intersection solver. Nonparallel nonaffine families, singular or collapsed intermediate sections, unproved angular or spatial embedding, unseparated loops and atlas components needing an unjoined self-boundary retain typed coverage refusals. The existing sketch restriction on mixing circular arcs with polynomial spline segments remains a sketch restriction; loft does not bypass profile validation. Elliptic sketch segments are outside this loft chart slice.

The two three-circle natural lofts above currently exhaust the eager application's unchanged 20M allowance. They remain explicit refusal regressions; this slice does not claim default-budget display or measurement for them. The next slice profiles their deterministic work before changing the implementation.

A separate qualification probe with circular radii `[1,2]` at heights `[0,3]` also returned `ResourceLimit` during `build_incremental_without_display`. No STEP or successful application result is claimed for that varying-radius circular family. It belongs with the work-cost follow-up, not with the successful analytic circle-prism cases.

Measurement enclosures retain their ordinary width marks. This feature does not change the numeric policy, modeling tolerance, generic joined-body Boolean classification or the separate rational-patch area-convergence follow-up. The field kernel does not gain loft evaluation from this B-rep entry.

SolidDraft still needs its Loft command and document conversion: collect profile references, one ordered start id per loop and section, and the ruled/smooth choice; emit the existing version 4/5 recipe fields. No SolidDraft source is changed here. Fixed spline sweeps are a separate operation slice.
