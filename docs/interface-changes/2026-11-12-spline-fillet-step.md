## 2026-11-12: spline fillet bodies in STEP

Pieces: shapes, STEP
Change: the STEP writer writes `SplineFillet` faces, faces on an atlas of one chart, and patch faces trimmed by a fillet's rails, with chart curves for every edge. `SplineFillet` gains `cut_curve`, which returns a `CutCurve`: a plane's cut as the carrier composed with a chart curve, with its proved bound. The writer adds the refusal code `step.spline-fillet-unproven`.
Why: M11 needs a rounded spline edge to open in OpenCASCADE ([prerequisite](2026-11-12-spline-fillet-body.md)).
Migration: none. Bodies the writer refused with `step.unsupported-geometry` for these faces now write or refuse with the new code.
