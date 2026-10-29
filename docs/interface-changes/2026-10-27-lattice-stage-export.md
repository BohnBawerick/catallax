## 2026-10-27: lattice stage export

Pieces: 5 (catallax-field `lattice_stage`, `printing`)
Change: adds `Stage::export`, `Output` and `Verdict`, as fixed in [2026-10-27-lattice-stage-contract.md](2026-10-27-lattice-stage-contract.md). `Settings.minimum_detail` now takes effect: boundary-coverage refinement may not use a grid edge below a quarter of it, and a request that needs one refuses as `field.print-unresolved` with a detail naming the minimum detail. The requested accuracy still sets the starting grid, so the proven bound is never loosened. `printing::mesh_plan` and every earlier refusal and byte are unchanged.
Why: lattice export must be closed, non-self-intersecting and connected at the user-chosen accuracy and detail.
Migration: none for existing callers. SolidDraft calls `Stage::export` and writes the mesh with `write_stl` or `write_obj`; the later immutable-certificate change below replaces direct mesh-field access with `Output::mesh()`.

`Verdict.non_self_intersecting` is true by construction (a closed 2-manifold on a conforming tetrahedral grid whose cells each hold one patch), not by a pairwise facet test at export time. The tests check it independently on a public synthetic lattice.

## 2026-10-30: labelled stage files and separate distance proof

Pieces: 5 (catallax-field `lattice_stage`, `printing`)

Change: adds `Output::labelled_stl(self) -> Result<(Vec<u8>, Output), Refusal>`. `Output` retains the proof's output controls privately. The method consumes the stage's proven mesh, writes the existing `catallax-print-stl/1` field header, and returns that mesh and its metadata-rounded verdict. The requested accuracy and caps come from the export, not from later caller edits. Header-range, encoding and `print.metadata-accuracy` refusals leave the destination untouched because the bytes are returned in memory.

Change: when accuracy is below the minimum-detail grid floor, stage export separates candidate meshing from two-sided distance certification. The candidate starts finer near bounding clip planes. All candidate and certificate cells, including failed refinements, share the supplied cap. The construction still proves closedness, non-self-intersection and one component, or refuses. See [FIELD.md](../FIELD.md#terminal-lattice-output) for the coverage argument. The general `printing::mesh_plan` route does not change.

Migration: obtain outputs through `Stage::export`; use `labelled_stl` for an active lattice's saved mesh. `Mesh::write_stl` remains the generic, unlabelled serializer and does not carry the stage's accuracy or caps.


## 2026-10-30: profile extrusion printing and immutable certificates

Pieces: 5, catallax-field `printing` and `lattice_stage`.

Change: printing accepts planar profile extrusions only after the existing exact convex-plane certificate proves their halfspace intersection. This covers lattices clipped to convex profile bosses without replacing their recipe with stock geometry. Other profile families retain `field.print-coverage`.

Change: `Output.mesh` and `Output.verdict` are private. Read them with `mesh()` and `verdict()`. Certified labels can no longer be written after mutating the certified geometry or verdict. The preview mesh remains mutable and unproved.

Migration: replace stage-output `.mesh` and `.verdict` reads with `.mesh()` and `.verdict()`. Clone the borrowed mesh if an ordinary, uncertified editable copy is needed. `labelled_stl` still consumes the original certified output.
