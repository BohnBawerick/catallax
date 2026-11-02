## 2026-11-02: all-faces references

Pieces: recipe language, face selection in both evaluators, SolidDraft conversion, verification driver.

Change: `catallax-recipe/5` permits an absent `FaceRef.face`. `{}` selects every face of the operation's final input. `touching` and `near` filter that full candidate set under their existing shared-edge and unique-nearest rules. A reference list remains a union. `[]` selects nothing; `[{}]` selects all. Named references remain required matches, not optional matches.

Migration: version 4 files migrate by explicitly changing `format` to `catallax-recipe/5`, checking and writing canonical text. Existing patterns and stored near points do not change. Version 5 inherits version 4 operation meanings. Review version 3 operations under the existing version 4 migration rules before upgrading. The parser keeps each saved version; versions 3 and 4 still reject references without `face`. `Recipe::new` retains its version 3 default so existing producers are not silently upgraded. Rust callers change `face: pattern` to `face: Some(pattern)`, or use `FaceRef::all()` for an unfiltered reference.

The legacy SolidDraft converter emits version 5 only when shell `open_faces: "all"` needs the new form. It writes `open: [{}]`, not references to every ancestor. This removes the structural-conversion refusal after booleans or finishing operations. Geometry coverage does not expand. A named `created_by` reference to a consumed producer still fails.

The independent through-bore field-shell construction also uses the optional origin filter. Unfiltered all-faces openings return the empty field, matching the B-rep result; selecting a bore wall for a nonempty shell retains its existing refusal.

### What SolidDraft must do

- Pin the merged Catallax change and accept/save version 5 in native recipe storage and requests. Keep older saved versions until an explicit migration.
- Write `[{}]` for an explicit all-faces choice, including shell openings, and `[]` for no selected faces. Do not enumerate face-producing ancestors or substitute the current body's creator origin.
- Keep picked faces as named references with `Some(pattern)` in Rust. Persist existing `touching` and `near` values unchanged through edits, undo, save and reopen.
- Retain typed geometry refusals. Do not treat all-faces syntax as new shell, draft or thicken geometry support. Unavailable controls show only `Coming soon`.
- Exercise all-faces selection after a boolean and a finishing operation, edit the source dimensions, save/reopen, undo/redo and verify the empty shell result and subsequent export handling. Also check that a vanished named producer still reports a lost selection. Catallax library checks do not establish rendered SolidDraft acceptance.

### Next slices, in order

1. SolidDraft native version-5 storage migration and all-faces picker acceptance.
2. Complete trimmed planar meeting arrangements and other shell openings.
3. Trimmed analytic and spline offset envelopes with matching shell and thicken domains.
