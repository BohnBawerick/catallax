## 2026-11-11: field preview reuses proved cell sides

Pieces: catallax-field display, catallax-app preview callers
Change: `catallax_field::display::preview_plan` no longer repeats the exact
witness checks on an open cell when classification has already proved which
side of the cell's support is material and that support belongs to a single
face. It takes the face and its orientation from that proof instead.
`PreviewMesh` gains `skipped_checks`, the number of face queries answered this
way. `mesh_plan`, `mesh_plan_reported` and every certified consumer are
unchanged and never skip a check.
Why: PLAN M13 F3. The certified request checks each emitting cell's witness
with exact rays along both normals and then evaluates the face's historical
ownership domain. Where the cell's two material sides are already proved and
no other face shares the support, those checks cannot change which face is
shown or how it is oriented. The [unproved hardware preview decision of
2026-11-10](../ARCHITECTURE.md#unproved-hardware-preview-evaluator-decision-of-2026-11-10)
keeps this result on the preview side: it is labelled unproved and never feeds
an export or a proven measurement.
Migration: code that builds a `PreviewMesh` by hand adds `skipped_checks`.
Supports shared by several faces, curved supports without a native fillet
chart, and cells whose sides are not proved keep the exact checks, so their
preview faces are chosen exactly as before. On P1 to P4 at 0.1 mm the preview
payload equals the certified payload; raw evidence retained privately.
