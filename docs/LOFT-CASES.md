# Rational loft cases

The rational loft entry is `catallax_app::build` with a version 4 or 5 `loft` operation. Its immutable definition is the recipe, including ordered start ids and `shape`. The correspondence law is [RECIPE 4.12](RECIPE.md#412-loft); this inventory follows [PLAN 3.8](PLAN.md).

Tests live in `crates/catallax-app/tests/rational_lofts.rs`. They construct actual recipes through the application entry, not manually assembled `Built` values.

The public successful recipes are [`loft-spline-prism.json`](../crates/catallax-features/tests/fixtures/loft-spline-prism.json) and [`loft-natural-circles.json`](../crates/catallax-features/tests/fixtures/loft-natural-circles.json). The second builds with display inside the unchanged 20M eager application allowance.

| Case | Contract test | Required result |
|---|---|---|
| Two-section ruled and smooth polygon lofts | `two_section_ruled_and_smooth_frustums_have_exact_mass_and_physical_topology` | Both laws give the same frustum; volume 28, centroid height 51/28, six physical faces, twelve edges and eight vertices. |
| Three circular sections, natural interpolation | `natural_three_circle_loft_builds_inside_the_default_allowance` | The saved recipe, radii 1, 2, 1 at heights 0, 1, 2, builds with display inside the unchanged 20M application allowance. Volume 383π/70, centroid height 1, three physical faces, two edges and no vertices; STEP face regions map to the model faces. |
| Unequal physical station distances | `unequally_spaced_natural_circles_build_inside_the_default_allowance` | The same radii at heights 0, 2, 5 build inside the same allowance. Uniform integer station spacing is retained, not chord-length interpolation: volume 383π/28, centroid height 43411/18384. |
| Two circles with different radii | `two_circles_with_different_radii_build_inside_the_default_allowance` | Radii 1 and 2 at heights 0 and 3 build inside the same allowance: volume 7π and centroid height 51/28, as for the exact frustum. |
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

A full circle's equal-angle quarter charts have algebraic weights `1, cos 45°, 1`. When every section's chart on a common parameter interval is a rational quadratic, Catallax applies one exact Möbius reparameterization to all of them: scaling the Bernstein weights by `λ^i` with `λ = w1/w2` keeps each chart's Euclidean controls and traces the same curve, with the new chart at `t` equal to the old chart at `λt/((1-t)+λt)`. That map fixes both interval ends, is strictly increasing and is the same for every section, so every correspondence pair `C_i(u)` stays together. The surface is the same point set with the same station curves; only the carrier's u speed changes. The weights become `1, k, k` with `k = cos²(span/2)`, which is 1/2 for quarter circles. The change applies only when exact field reduction proves every new coordinate rational and every new weight positive; otherwise the interval keeps its algebraic chart. A loop whose charts are all rational then carries a running rational homogeneous scale across its intervals instead of raising the u degree. The unit test `one_common_mobius_map_makes_quarter_circles_rational_without_moving_them` checks the map exactly at interior parameters.

Exact coefficient identities reduce straight translated section families to existing prisms when their profile carriers permit it. Positive rational control-hull proofs reduce planar logical sides and straight physical boundaries to ordinary planes and lines. These are exact reductions of the recipe definition, not fits, and avoid running rational-surface integration and display on analytic geometry.

These certificates are sufficient conditions, not a general loft intersection solver. Nonparallel nonaffine families, singular or collapsed intermediate sections, unproved angular or spatial embedding, unseparated loops and atlas components needing an unjoined self-boundary retain typed coverage refusals. The existing sketch restriction on mixing circular arcs with polynomial spline segments remains a sketch restriction; loft does not bypass profile validation. Elliptic sketch segments are outside this loft chart slice.

Measurement enclosures retain their ordinary width marks. This feature does not change the numeric policy, modeling tolerance, generic joined-body Boolean classification or the separate rational-patch area-convergence follow-up. The field kernel does not gain loft evaluation from this B-rep entry.

SolidDraft still needs its Loft command and document conversion: collect profile references, one ordered start id per loop and section, and the ruled/smooth choice; emit the existing version 4/5 recipe fields. No SolidDraft source is changed here. Fixed spline sweeps are a separate operation slice.

## Where the 20M allowance goes

This count follows the deterministic work profile of [NOT-EXACT-INVENTORY.md](NOT-EXACT-INVENTORY.md), applied to the eager application build. The `app_work_profile` example (feature `work-profile`) runs `catallax_app::build` on one recipe with display, chord tolerance 1e-2 and 200,000 display vertices, the suite's options. It attributes every unit of the build's work to a decider call, its calling line and the logged operation. Work outside the loft operation is display and naming. Work units are the same on any machine and under any load.

```text
cargo run --release -p catallax-cli --features work-profile --example app_work_profile -- RECIPE --max-work UNITS
```

The three recipes are the suite cases above. Raw evidence is retained privately.

### Before: algebraic quarter-circle charts

| Case | 20M result | Larger allowance |
|---|---|---|
| Two circles, radii 1 and 2 | `ResourceLimit` inside the loft operation | Builds at 4G: 182.0M in total, 46.5M in the loft operation and 135.5M in display |
| Three circles, heights 0, 1, 2 | `ResourceLimit` inside the loft operation | Still `ResourceLimit` inside the loft operation at 4G |
| Three circles, heights 0, 2, 5 | `ResourceLimit` inside the loft operation | Still `ResourceLimit` inside the loft operation at 4G |

At 20M the three-circle lofts spend 11.7M enclosing the projection frame of the rational-patch embedding proof (`catallax-shapes/src/math.rs`), 7.0M enclosing B-spline controls (`catallax-shapes/src/bspline.rs`) and 1.0M evaluating polynomials with algebraic coefficients. Each control is an expression over `√(1/2)`, and every enclosure walks its whole definition. At 4G, 3.93G is exact-zero proofs that adjacent chart cells share their seam (`RationalPatch::seam`, `catallax-shapes/src/rational_patch.rs`): 672 calls at about 5.8M units each, and the operation is still unfinished. The two-circle loft finishes at 4G. Its loft operation needs 46.5M, of which 26.2M is the same seam proofs; its display needs 135.5M, including 39.4M of exact-zero proofs in planar display (`catallax-mesh/src/planar.rs`) and 39.4M of edge enclosures (`catallax-mesh/src/edges.rs`). No small saving at one call site closes that gap: the cost is in the algebraic carrier.

### After: one common rational reparameterization

| Case | Total at 20M | Loft operation | Display and names | Display triangles |
|---|---|---|---|---|
| Two circles, radii 1 and 2 | 9,174,581 | 522,857 | 8,651,722 | 38,248 |
| Three circles, heights 0, 1, 2 | 17,401,490 | 1,328,542 | 16,072,945 | 113,180 |
| Three circles, heights 0, 2, 5 | 18,551,052 | 1,228,198 | 17,322,851 | 125,468 |

The two sketches or three sketches add one unit each. In the three-circle cases, B-spline evaluation in display (`catallax-shapes/src/bspline.rs`) takes 14.5M to 15.7M: 6.4M to 6.7M of control enclosures and 8.1M to 9.0M of deadline checks, one unit each. The rational-patch display grid takes about 0.6M. In the two-circle case, B-spline evaluation takes 4.1M and sign filters in `catallax-shapes/src/rational_patch.rs` 1.7M.

The loft operation now costs between 0.5M and 1.4M units. Display is the remaining cost, and each control enclosure there is a fixed charge on a rational value. The three-circle cases keep a margin of 1.4M to 2.6M units. No allowance, tolerance, display cap or work charge changed. The 20M allowance is still the application default, and each case reports `ResourceLimit` if it is exhausted.
