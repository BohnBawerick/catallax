## 2026-10-30: export-only closing of a planar wedge too narrow to read

Pieces: 7, CLI verification

Change: a planar face whose wire turns back across one straight edge with a
written length at or below 1e-7 mm, between two straight edges at least ten times
longer that border planar walls facing the same way, is closed in the file. One
wall moves onto the other's plane within the requested tolerance and 1 nm, the
short edge and both long edges leave the file, and one straight edge joins their
far ends unless they are one vertex. Each record reports one removed face and two
removed edges, or two faces and three edges for a three-sided wedge face. A record
that follows an earlier repair can name a face or edge of the writer's copy rather
than a model handle. `StepRegularization` and `catallax-step-regularizations/1`
keep their fields and meanings, and follow
[the narrow-step collapse](2026-10-30-step-narrow-step-collapse.md).

Why: SolidDraft's GUI frames lean the bracket's plate front 2e-17 mm proud of the
flange and start the flange 2e-28 mm in from the plate's side. The exact union
correctly keeps a wedge in front of the flange that closes to a 2e-17 mm edge, so
the file was raw-invalid from the first union onward, and the raw reader returned
two zero-area faces.

Migration: the exact model, names, marks, measured properties and manifest counts
do not change. Consumers that add every record's `removed_faces` to the file's
face count need no change. Files with no such wedge keep their bytes.
