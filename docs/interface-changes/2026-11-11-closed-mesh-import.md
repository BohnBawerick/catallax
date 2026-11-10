# 2026-11-11: closed mesh import in both kernels

Pieces: import, features, field, compare, application

PLAN M22. `import_mesh` (RECIPE section 4.19) now builds in the B-rep and in the field. No existing outcome, default or saved answer changes; recipes without `import_mesh` take the same paths as before.

## New interface

- `catallax_import::mesh::ClosedMesh` reads binary STL, ASCII STL or Wavefront OBJ bytes with an exact positive scale and proves, with no tolerance, that the triangles bound a solid. It exposes the exact vertices, triangles, edge twins, shells with orientation and nesting, lumps, exact bounds, the integer winding number of a point, and the exact squared distance from a point. `catallax_import::mesh::MAX_TRIANGLES` is the stated limit, 200,000 triangles.
- `catallax_features::ImportFiles` now also carries mesh bytes, keyed by the SHA-256 saved in each `import_mesh` step. `evaluate_with_files`, `evaluate_incremental_with_files` and `catallax_app::build_incremental_with_files` build the step. Missing or changed bytes are recipe errors, worded as for STEP sources.
- `catallax_field::Plan::compile_with_files(recipe, files, dec)` compiles with the same byte map. `Plan::compile_checked` is that call with no files, so it refuses an `import_mesh` step with `field.recipe`.

## Refusal codes

Impossible: `import.mesh-open`, `import.mesh-orientation`, `import.mesh-degenerate`, `import.mesh-crossing`, `import.mesh-winding`, `import.mesh-scale`. Unsupported: `import.mesh-nonmanifold-edge`, `import.mesh-pinched-vertex`, `import.mesh-too-large`. Unreadable bytes: `import.mesh-malformed`. Both kernels return the reader's code unchanged; the B-rep wraps it in `BuildError::Geometry` for the import step.

## Results

The B-rep body has one planar face per edge-connected set of coplanar triangles, every face named `{step, role: "mesh"}`, one shell per connected surface and one lump per outer shell with its cavities. Its exact volume equals the field's. The field leaf is true signed distance; a lone mesh has exact field volume, centroid and inertia. The field registers no named mesh face, a ray from a boundary point stays `Indeterminate`, and field mass of a Boolean holding a mesh refuses with `field.print-coverage`. The comparison harness reports `import_mesh` as not comparable because it carries no source bytes.

## Migration

Callers that build recipes with `import_mesh` pass the mesh bytes in the existing file map. Field callers switch from `Plan::compile_checked` to `Plan::compile_with_files` for such recipes. Nothing else changes.

## Case inventory

Tests are in `crates/catallax-import/tests/mesh.rs`, `crates/catallax-features/tests/import_mesh.rs` and `crates/catallax-field/tests/import_mesh.rs`.

| Case | Test |
| --- | --- |
| Closed cube, winding through edges and vertices, a 10^-30 offset, exact distance | `a_closed_cube_bounds_its_interior_by_winding_number` |
| The same exact mesh from binary STL, ASCII STL and OBJ | `binary_ascii_and_obj_files_give_the_same_exact_mesh` |
| Exact binary32 coordinates, exact scale, nonpositive scale | `binary_coordinates_are_the_exact_binary32_values`, `scale_multiplies_exact_file_coordinates` |
| Open mesh refused with the open edge | `an_open_mesh_is_impossible_and_names_the_open_edge`, `open_and_crossing_meshes_are_refused_as_impossible_with_the_reason`, `the_field_refuses_open_crossing_and_missing_meshes` |
| Reversed triangle, inside-out surface, cavity facing the material | `a_reversed_triangle_is_an_orientation_error`, `an_inside_out_surface_has_winding_minus_one_and_is_refused`, `a_cavity_must_face_into_the_void` |
| Crossing bodies, touching faces, a folded pair, a corner pushed through a face | `overlapping_bodies_cross_and_are_refused`, `bodies_touching_along_a_face_cross`, `a_folded_pair_of_triangles_crosses`, `a_surface_folded_through_itself_crosses` |
| Zero-area triangles; shared edges and pinched vertices as unsupported | `zero_area_triangles_have_no_orientation`, `shared_edges_and_pinched_vertices_are_unsupported_not_impossible` |
| Cavities, islands and separate bodies as shells and lumps | `separate_bodies_and_an_island_inside_a_cavity_are_separate_lumps`, `a_cavity_is_an_inner_shell_and_a_hole_is_an_inner_loop` |
| Coplanar neighbours merged into planar faces, hole loops, genus one | `a_finely_split_flat_face_is_valid`, `a_closed_mesh_builds_an_exact_solid_with_merged_planar_faces`, `a_cavity_is_an_inner_shell_and_a_hole_is_an_inner_loop` |
| Scale and pose in the B-rep | `scale_and_pose_place_the_exact_solid` |
| Field classification and true distance | `the_field_decides_inside_by_winding_number_with_true_distance` |
| Exact volume, centroid and inertia equal in both kernels | `both_kernels_give_the_exact_polyhedron_volume`, `scaled_rotated_meshes_agree_exactly_between_kernels`, `a_curved_binary_mesh_has_one_exact_volume_in_both_kernels` |
| Size limit with a typed refusal in the reader and both kernels | `the_size_limit_is_a_typed_refusal`, `a_mesh_above_the_stated_limit_is_refused_in_both_kernels` |
| Missing or changed bytes | `missing_or_changed_bytes_are_recipe_errors_and_zero_scale_is_impossible`, `the_field_refuses_open_crossing_and_missing_meshes` |
| Malformed files and the empty mesh | `malformed_files_say_why`, `an_empty_closed_mesh_bounds_nothing`, `an_empty_mesh_builds_the_empty_solid`, `an_empty_mesh_is_the_empty_field` |
| A mesh in a field Boolean | `a_mesh_in_a_field_boolean_classifies_but_has_no_mass_coverage` |
