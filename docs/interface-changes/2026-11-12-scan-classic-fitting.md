## 2026-11-12: Classic fitting of scans to a clean part

Pieces: field, mesh and scan import

Change: `catallax_field::reconstruction::scan::fit_scan` fits a `scan_import` step's point cloud or mesh with region growing and least-squares plane, cylinder and cone candidates, builds an ordinary `catallax-recipe/4` recipe (the scan step, then a block stock, bosses, bores and stock-edge rounds) in both kernels, and reports both sampled deviation numbers with their maps. New refusal codes: `reconstruction.unexplained-region`, `reconstruction.unplaced-surface`, `reconstruction.no-stock`, `reconstruction.too-few-samples`, `reconstruction.too-many-samples`, `reconstruction.recipe`, `reconstruction.invalid-part`, `reconstruction.nonfinite-deviation`, `reconstruction.nonfinite-sample` and `reconstruction.deviation-too-large`. `catallax_import::scan::FastSurface` exposes the unproved binary64 nearest-distance index the scan reference already used.

Resource checks: tree construction, visited search nodes and candidates, fitting iterations, graph visits and deviation statistics charge the caller and check its deadline. Sample counts are checked before allocation, with at most 1,000,000 samples per deviation direction. `FastSurface::new_checked` and `distance_checked` expose the metered index path; existing unchecked methods retain their behavior.

Why: PLAN M23, classic fitting comes first. Details, thresholds and limits are in [M23-SCAN-FIT.md](../M23-SCAN-FIT.md).

Migration: none. The operation is opt-in; no existing route calls it, and no saved answer, tolerance or allowance changes.
