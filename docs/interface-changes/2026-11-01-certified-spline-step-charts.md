## 2026-11-01: Certified spline prism STEP side becomes two exact charts

Pieces: 7, judge interface
Change: `write_spline_prism_step` writes the one model side as two open B-spline faces on `[0,2]x[0,1]` and `[2,4]x[0,1]`, with two half arcs per ring and two shared axial edges. The file has 4 faces, 6 edges and 4 vertices. `StepReport` gains `face_map: Vec<StepFaceRegion>`, empty for other writers. `step.regions` of a nonempty certified spline answer has four rows, and the two side rows add an exact `file_domain`. `step.added_vertices` goes from 2 to 4 and `step.added_edges` from 0 to 4.
Why: the closed side had two ring edges and no seam, so direct no-healing translation gave an invalid solid and a zero-area side. The split keeps the polynomial geometry exactly: bounds, certificates and native counts do not change.
Migration: read a file face's name from `step.regions` or `StepReport.face_map`, not from a model handle. Two rows can share one `face`. Regenerated answers from this writer have new bytes. Old raw-invalid files stay as they are.
