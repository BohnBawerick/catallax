## 2026-10-30: export-only removal of a planar step too narrow to read

Pieces: 7, CLI verification

Change: `StepRegularization` adds `removed_faces` and `removed_edges`, and
`step-regularizations.json` adds the same two keys to each record. The format name
stays `catallax-step-regularizations/1`, because the keys are additive. A four-edge
planar ribbon with written end edges at or below 1e-7 mm is removed from the file:
one neighbouring wall moves onto the other's plane within the requested tolerance
and 1 nm, and the two walls become one face. Each such repair removes two faces and
four edges from the file. The record for an existing wire-gap reconnection reports
one face and no net edges.

Why: OpenCASCADE cannot translate an edge at or below 1e-7 mm, whatever the file
declares. SolidDraft's plate and flange at width 120 meet with a real 2.8e-15 mm
step, and the exact body correctly keeps it. The file was raw-invalid, and from a
width near 5e-9 mm up to 1e-7 mm the standard reader hid that by deleting the short
edges and losing volume.

Migration: the exact model, names, marks, measured properties and manifest counts
do not change. A consumer that compares file and model face counts must add every
record's `removed_faces` to the file's count. SolidDraft's `ci/check-step.py` needs
that change. Files with no such ribbon keep their bytes.
