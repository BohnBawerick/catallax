## 2026-11-09: Checked finite rational surface patches

Pieces: 2, 3, 4, 5, 6, 7; display mesh and application integration tests.

Change: `catallax_shapes::RationalPatch` adds a checked exact finite rational surface carrier, represented by `SplineSupport::RationalPatch`. Complete rectangular faces have finite predicates, display meshes, measurement and bounded STEP export. This does not implement loft/sweep feature builders or expose new recipe/application commands.

Why: generic rational carriers and their consumers must land before stage 5 operation builders emit them. The [versioned path/correspondence schema](2026-11-10-rational-path-schema.md) is a separate prerequisite; this change does not modify that schema.

Migration: callers constructing this support must first validate a `BSplineSurface`, then call `RationalPatch::new`. Use `RationalPatchFace::of_face` before consuming a complete physical face. `catallax_intersect::supports` now takes `&Surface` arguments rather than `SurfaceKind` arguments: a coarse spline kind cannot authorize a particular spline family. Its in-tree caller is migrated; the old signature is not retained.

### Exact definition and proof boundary

The clamped knots, controls and positive weights remain authoritative. Construction changes none of them. Exact span equations support queries and integration. `BSplineSurface::span_bezier_controls` and `BSplineCurve::span_bezier_controls` share convex blossom extraction: nonuniform patch rims are bounded as actual Bézier spans, not raw B-spline control polygons. No sampled surface replaces the definition.

A complete-cell positive symmetric-Jacobian bound establishes a globally injective regular projection. The proof uses unit-domain coordinates, independent of the affine units of authored knots. `domain`, `point_at`, `derivatives_at` and shared-row stations retain original knot coordinates. `embedding_lower_bound` and `embedding_cover` expose the certified unit-domain projection margins. Internal position-continuous knot joins must also have compatible oriented tangent planes. A failed certificate returns `shapes.rational-patch-unproven`, not an impossible-geometry claim.

`RationalPatchFace` checks four complete physical rows, endpoint identity, a sewn vertex cycle and material-relative orientation. It does not certify nonintersection between different faces of a body. Arbitrary trimming, holes in the parameter rectangle and intersection-defined boundary correspondence remain outside this face certificate. Matching endpoints alone never authorizes a row. Names stay on the existing physical faces, edges and vertices; reflection preserves origins and reverses the tensor face/loop orientation as required.

Finite preimages validate original world-coordinate equations. Normal and side queries do not extend the finite chart into an infinite implicit surface. Whole-line intersections enumerate transverse finite hits; tangencies and contained curve families refuse. Plane sections retain their generating patch and plane, and emit exact rational curves for the supported univariate/affine-parameter cases. General unresolved algebraic tracing remains a typed gap.

### Mesh, measurement and file consumers

Display grids include knots and exact required shared stations. Adjacent rows exchange physical points through exact inverse parameters, including nonlinear reparameterisation; they do not assume matching station indices. Complete-cell Hessian bounds certify barycentric chord error. Unit-normal cones remove parameter-speed variation using `|d(n/|n|)| = |n cross dn| / |n|²`. The mesh writer separately checks its rounded facets against those cones and verifies rounding against the sealed projected embedding margin. Display error is not a geometry-fitting tolerance.

Measurement integrates unit-local exact Bézier definitions, with validated tensor Gauss remainders in both parameter directions. The homogeneous normal cancels its known positive weight factor symbolically before interval evaluation. Volume, area, centroid, inertia and extrema use the retained definition, never the display mesh. The existing accuracy-mark contract remains: an exhausted request keeps its certified enclosure and reports the unmet width.

STEP export writes rational surfaces, shared three-dimensional rows and their finite pcurves. Authored knot coordinates normalize affinely for the file; decimal control, weight and knot errors have a separate reported file bound. Shared rows require a proved common normalized chart. A nonlinear chart mismatch returns `step.rational-patch-unproven` without publishing partial output, even when display meshing can share the physical row. Reading exported source data does not add certified generic-spline solid import support.

### Known gaps and follow-ups

Area integration retains truthful marked enclosures when the requested width cannot be reached under the default work budget. No numerical budget or geometry tolerance is changed to hide this limit. Separate mass tests request reachable area accuracy while retaining independent closed-form area, volume, centroid, inertia and body-bound truth checks. Follow-up: improve area convergence within the existing budget.

| Area case | Requested absolute width (mm²) | Observed enclosure width (mm²) | Outcome |
| --- | ---: | ---: | --- |
| Bilinear roof solid | 1e-7 | approximately 1.604629e-7 | `measure.accuracy-not-reached`, with the enclosure width as its value bound |
| Nonuniform quadratic bump solid | 1e-6 | approximately 4.202163e-3 | `measure.accuracy-not-reached` |

Joined degree-elevated and joined nonlinear fixtures retain `boolean.parametric-query` for body classification: the relevant Boolean path still needs a finite parametric curve solver. Individual complete-face witnesses/membership and shared-edge display work. Follow-up: close that finite Boolean gap without treating a finite spline arc as an unbounded line. This change does not add a line-closure reduction or another Boolean solver.

Reachable-accuracy checks request absolute measurement accuracy 1e-6 for the bilinear/rational solids and 1e-2 for the nonuniform bump. Volume, centroid and inertia enclosures still must meet the independent 1e-6 width checks, and every variant's body bounds are requested and checked at 1e-6. Both area definitions are compared with certified closed-form oracles, not sampled reference values.

### Case inventory

Shape cases live in [`catallax-shapes/tests/rational_patches.rs`](../../crates/catallax-shapes/tests/rational_patches.rs). Consumer cases and their closed named fixtures live in [`catallax-app/tests/rational_patches.rs`](../../crates/catallax-app/tests/rational_patches.rs) and its `rational_patch_cases` module. These are public-kernel tests, not claims that a loft/sweep recipe reaches application display.

| Case | Test / required outcome |
| --- | --- |
| Positive rational weights; exact point and finite inverse | `positive_rational_weights_keep_exact_positions_and_finite_incidence` |
| Near-but-off-surface point; finite-domain exclusion | Same test: proven absence, not snapping |
| Transverse line hit and proven finite miss | `transverse_line_roots_are_complete_and_tangency_is_not_a_false_miss` |
| Tangent line and contained line family | Same test: typed refusal, not an empty hit list |
| Rational hyperbolic plane section | `hyperbolic_plane_sections_retain_the_exact_curve_and_definition` |
| Folded definition and internal crease | `folded_and_creased_definitions_cannot_acquire_an_embedding_certificate` |
| Flat G1 parameter-speed jump, including angular control | `a_g1_speed_jump_keeps_angular_mesh_control` |
| Huge translated, non-dyadic knots | `huge_translated_non_dyadic_knots_do_not_lose_the_chart` |
| Interior extrema distinct from the control polygon | `exact_ranges_enclose_geometry_rather_than_the_control_polygon` |
| Complete bilinear/rational/bump faces; exact witnesses, membership and off-surface exclusion | `complete_named_faces_have_exact_witnesses_and_finite_membership` |
| Joined degree-elevated and nonlinear body classification | Same test: explicit `boolean.parametric-query` gap |
| Tight bilinear area request | `tight_bilinear_area_keeps_truth_and_reports_the_unmet_width`: contains the independently enclosed closed-form area; accurate value mark |
| Bilinear/rational/bump mass and complete nonuniform-rim body bounds at reachable requested accuracy | `rational_parameterisation_preserves_mass_moments_and_geometric_bounds` |
| Degree changes and nonlinear shared-row parameters | `shared_rows_sew_even_with_degree_changes_and_nonlinear_parameters`: closed oriented mesh, physical face keys and certified controls |
| Normalized rational STEP surfaces and shared charts | `step_rounding_preserves_the_source_equations_and_shared_row_charts` |
| Nonlinear shared STEP chart | Same test: typed refusal and no partial bytes |
| Warm definition caches followed by reflection | `warmed_definition_caches_do_not_leak_through_reflection`: names, positive volume and transformed centroid |
| Nonzero tensor Gauss remainder in both directions | `integrate::gauss::tensor::tests::product_rule_bounds_a_nonzero_remainder_in_both_directions` |

Native STEP qualification uses generated kernel fixtures, an OpenCASCADE reader with explicitly disabled shape processing after `ReadFile`, raw topology checks and independent analytic mass properties. Pinned replay compatibility remains a byte-for-byte gate; neither known gap authorizes changing replay answers.

SolidDraft is unchanged. Its later Loft/Sweep commands need operation builders, application dispatch, version-4 wire conversion, explicit loft correspondence and the documented fixed-path frame definition before presenting these carriers as new modeling commands.
