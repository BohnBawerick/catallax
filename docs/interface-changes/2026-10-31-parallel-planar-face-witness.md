## 2026-10-31: parallel planar face distance with exact gap witnesses

Pieces: 6 (measurement), SolidDraft Catallax adapter

Change: Added `PreparedFace::parallel_face_to_face`. For planar faces with exact rational frames and parallel normals, the supporting-plane gap is the lower bound, and an exact probe admitted by both trims attains it. Without such a witness it returns the general `face_to_face` bounds. Curved faces, nonparallel planes and inexact frames refuse `measure.nonplanar-face`, `measure.nonparallel-faces` and `measure.irrational-plane`. `face_to_face` is unchanged.

Why: At 1e-8 mm the general search stopped at width about 3.3e-5 mm after 4,095 cells on the holed caps after a radius 2 fillet of one bottom edge. The witness path answers those caps exactly.

Migration: Existing callers need no changes. Face pairs within the parallel planar domain can call the new method with an explicit `Options::width` and must still check `converged` and `upper`.
