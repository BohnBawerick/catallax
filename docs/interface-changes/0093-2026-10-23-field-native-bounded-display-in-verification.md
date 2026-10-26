## 2026-10-23: field-native bounded display in verification

Pieces: field, verification CLI, verification driver
Change: `catallax_field::display` exports `Options`, `Mesh`, `Triangle` and `mesh_plan`. Field requests construct their own certified display mesh and write `display.obj`, with `face_measurements.display.first_mesh_ms` measured across recipe processing, field compilation and native meshing. The display report adds a two-sided boundary-distance upper bound, artifact name and deterministic cell count. The driver accepts `--source field` and compares strict OBJ bytes as well as STEP.
Why: plan metadata alone is not a first display, and B-rep construction cannot measure a field-first feature.
Migration: field callers must handle display refusals as well as compile refusals. Display limits default to 0.01 mm, one million vertices and two million parameter cells. New refusals are `field.display-options`, `field.display-unsupported`, `field.display-unbounded`, `field.display-budget` and `field.display-precision`. Refusals remove stale OBJ and STEP artifacts. Unsupported operations retain `field.unsupported-operation`; no B-rep fallback is added. Display bounds are not whole-face fillet construction tolerances. Legacy and B-rep answers, recipe construction and STEP options are unchanged.

