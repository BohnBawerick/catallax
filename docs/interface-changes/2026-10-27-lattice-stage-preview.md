## 2026-10-27: unproved lattice preview

Pieces: 5 (catallax-field `lattice_stage`), judge interface (`field` verification adapter)
Change: adds `Stage::preview(&PreviewOptions, dec) -> Result<Preview, Refusal>`, `Preview { mesh, proof: Proof::Unproved }`, `lattice_stage::preview_plan` and `lattice_stage::default_preview_options`. This replaces the `Stage::preview(&self, dec)` signature in [2026-10-27-lattice-stage-contract.md](2026-10-27-lattice-stage-contract.md). The preview is the print-grid mesher run at the preview's own tolerance and limits, from the field alone. It is not the export, carries the coarse grid's own `max_deviation`, and requires neither connectivity nor the user's accuracy.
Change: the `field` verification adapter, when the certified display refuses with `field.display-lattice`, now answers with that unproved preview: `exact` is `false` and `face_measurements.display.proof` is `"unproved"`. Every other field request, its limits and its refusals, including all other uses of `field.display-lattice` outside the adapter, are unchanged. The preview's refusals (`field.print-budget`, `field.print-grid`, `field.print-vertices`) are returned as before.
Why: the lattice stage needs a visible display before proof.
Migration: callers of the field adapter that expected the `field.display-lattice` refusal for a lattice recipe now receive an unproved answer and must not treat it as exact.

## 2026-10-30: sampled cubic preview

Constant axis-aligned cubic lattices now sample the existing tentative field on a bounded rectilinear grid rather than running the certified print mesher. The pitch uses the preview target and beam radius, independently of export settings. Cell and vertex exhaustion still refuse without returning a partial mesh. Other lattice previews retain the previous interval-driven route.

The result remains `Proof::Unproved`. Its mesh's `max_deviation` reports a sampling scale only, not a proven distance bound. Callers must keep the unproved label and obtain proof separately through `Stage::export`.

