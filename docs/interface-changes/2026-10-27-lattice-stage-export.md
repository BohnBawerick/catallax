## 2026-10-27: lattice stage export

Pieces: 5 (catallax-field `lattice_stage`, `printing`)
Change: adds `Stage::export`, `Output` and `Verdict`, as fixed in [2026-10-27-lattice-stage-contract.md](2026-10-27-lattice-stage-contract.md). `Settings.minimum_detail` now takes effect: boundary-coverage refinement may not use a grid edge below a quarter of it, and a request that needs one refuses as `field.print-unresolved` with a detail naming the minimum detail. The requested accuracy still sets the starting grid, so the proven bound is never loosened. `printing::mesh_plan` and every earlier refusal and byte are unchanged.
Why: lattice export must be closed, non-self-intersecting and connected at the user-chosen accuracy and detail.
Migration: none for existing callers. SolidDraft calls `Stage::export` and writes `Output.mesh` with `write_stl` or `write_obj`.

`Verdict.non_self_intersecting` is true by construction (a closed 2-manifold on a conforming tetrahedral grid whose cells each hold one patch), not by a pairwise facet test at export time. The tests check it independently on a public synthetic lattice.
