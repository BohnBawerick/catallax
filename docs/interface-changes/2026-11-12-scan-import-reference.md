## 2026-11-12: `scan_import` reads open scans as reference data

Pieces: recipe format, B-rep recipe adapter, field recipe adapter, comparison harness, STEP and mesh import

Change: `catallax-recipe` gains `Operation::ScanImport` (`scan_import`, version 4 and later) and `ValueKind::Scan`. `catallax_import::scan::ScanReference` reads the step's bytes with the mesh readers, keeps holes, flipped and crossing triangles, degenerate triangles and isolated samples, places them exactly and reports what the file holds. New refusal codes: `import.scan-too-large`, `import.scan-scale`, `import.scan-hash` and `import.scan-step`. The mesh reader keeps unused OBJ vertices as isolated samples for scans; closed mesh import ignores them as before.

Why: PLAN M23 starts from scan import as reference data, holes and all, as RECIPE.md section 4.19.1 defines it.

Migration: an exhaustive `match` on `Operation` or `ValueKind` needs the new arm. A scan makes no body: the B-rep and field adapters skip the step, and the comparison harness has nothing to pair. The recipe check refuses a scan as the result, as a body input and with `scale <= 0`. No existing recipe, saved answer, tolerance or allowance changes.
