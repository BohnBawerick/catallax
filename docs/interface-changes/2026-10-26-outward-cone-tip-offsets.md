## 2026-10-26: outward cone-tip offset sectors and algebraic cone witnesses

Pieces: features, field, display mesh

Change: outward perpendicular-capped cone offsets directly assemble the spherical normal sector at the apex and its shared native contact circle on the parallel cone support. `ConeSkin::outward_pieces` unions the finite pointed-wall normal image, the tip sphere and the base edge torus independently of B-rep topology. `ConeSweep::new` reduces rational meridian differences and axial apex shifts before building cone supports. Field distance queries and candidate charts evaluate algebraic parallel cones at the rational apex section and along rational stereographic generators without introducing nested radial square roots. `FaceBuilder::finish` scales a cone's second-derivative $u^2$ term by each triangle's own maximum axial parameter $|v|$ rather than the whole-face maximum, matching `Chart::edge_chord_within`.

Why: pointed outward cone offsets need a spherical normal sector at the apex, radical-free field witnesses when the cone slope is not a Pythagorean ratio, and a triangle-local cone Hessian bound when the base-to-tip rail radius ratio exceeds eight.

Migration: none. Existing recipe formats, tolerances, work allowances and earlier answers are unchanged. Pointed conical band sweeps, oblique outward conical edges, coaxial cylindrical neighbours and non-ring edge tubes retain typed refusals.
