## 2026-10-27: final lattice stage contract

Pieces: 5 (catallax-field `lattice_stage`)
Change: adds `catallax_field::lattice_stage` with `LatticeSpec`, `Settings`, `Stage` and `check_no_cad_after_lattice`. Merged with this note: request building, rebuild, return before the lattice, and the CAD-after-lattice refusal. Preview, export and the verdict arrive in later slices with the signatures below.
Why: SolidDraft makes the lattice the last feature of a finished body and needs one stable call surface for it.
Migration: none for existing callers. SolidDraft builds against the section below.

### Merged now

`Stage::new(body: &Recipe, spec: LatticeSpec, settings: Settings, dec: &mut Decider) -> Result<Stage, Refusal>`

- `body` is the finished body as a `catallax-recipe` (SolidDraft converts its history first). Its result step is the body the lattice fills. The stage appends one lattice step reading that result and makes it the result; `Stage::recipe()` and `Stage::lattice_step()` return them. The step id is `lattice`, or `lattice-N` if the body already uses that name.
- `LatticeSpec { cell, size, frame, radius }` is the lattice the user chose. `Settings { accuracy, minimum_detail, max_cells, max_vertices }` are the output controls in millimetres. Only accuracy and minimum detail are meant to be relaxed by the user; the stage never changes the body, its tolerances or the 1 nm modelling tolerance.
- A `Stage` holds no cache and no state shared with an earlier body. Editing the body means calling `Stage::new` again with the edited recipe, the saved `LatticeSpec` and the saved `Settings`. Dropping the old `Stage`, and any preview or output made from it, is the whole invalidation.
- `Stage::return_before_lattice(self) -> (Recipe, LatticeSpec, Settings)` gives back the unchanged body recipe with the lattice parameters and settings, so SolidDraft can keep them while the user edits.
- `check_no_cad_after_lattice(&Recipe) -> Result<(), Refusal>` reports a history in which any step reads a lattice.

Typed refusals (`Refusal.code`):

| Code | Meaning |
|---|---|
| `lattice.cad-after-lattice` | The body already contains a lattice, or a step reads one. Return before the lattice. |
| `lattice.settings` | Accuracy, minimum detail or a limit is not positive, or fewer than four vertices are allowed. |
| `lattice.impossible` | The cell size, cell or beam radius is impossible. |
| `field.*` | The body cannot be compiled by the field; the field's own code and detail are passed through. |

### Next slices (signatures fixed here)

- `Stage::preview(&self, dec) -> Result<Preview, Refusal>`: an unproved display of the lattice, reported as `Proof::Unproved`. It never uses B-rep geometry. Typed refusals are the field's own (`field.display-*`, `field.print-budget`).
- `Stage::export(&self, dec) -> Result<Output, Refusal>`: the certified mesh at `Settings.accuracy`, with a `Verdict { closed, non_self_intersecting, connected, components, max_deviation }`. Output is only returned when all three are true. Refusals: `field.print-disconnected` (disconnected), `field.print-unresolved` and `field.print-topology` (not provable), `field.print-budget`, `field.print-grid` and `field.print-vertices` (over budget). Deadline, I/O failure and cancellation stay infrastructure outcomes, not geometry.
- `Settings.minimum_detail` will bound refinement: where the requested accuracy needs a finer grid than this detail allows, the export refuses as not provable instead of dropping the feature.
