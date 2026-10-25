## 2026-10-25: seed fine revolution meshes from certified shared rails

Pieces: display mesh

Change: eligible cone and bounded-torus interiors use the measured maximum longitude gap in their shared boundary samples. Only the remaining parameter subdivides when that fixed longitude already has a strict chord and angular margin. Sparse trims retain the general grid. Interior columns merge chart-rounding duplicates without modifying shared samples; apex rows retain the existing canonical apex handling. Every returned facet still passes the complete chord and normal certificates.

Why: a balanced cone grid ignores the much finer longitude already present on its boundary, and a bounded torus's nominal grid spacing can overestimate the actual shared-rail gap. Both generate avoidable interior samples and triangulation work.

Migration: none. Requested tolerances, user resource limits, exact definitions, STEP output and shared-edge ownership are unchanged. Deterministic display connectivity and vertex counts can change. Rotated-cone regressions compare compact output counts with the measured d0cea699 baseline at identical controls, rather than replacing the temporary-chart budget with an arbitrary output limit.

Prior art: [OCCT's cone range splitter](https://github.com/Open-Cascade-SAS/OCCT/blob/master/src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_ConeRangeSplitter.cxx) scales parameter steps using cone radius and excludes boundary rows. [Shewchuk's anisotropic interpolation discussion](https://people.eecs.berkeley.edu/~jrs/papers/elemj.pdf) separates interpolation accuracy from isotropic element shape. This change keeps Catallax's own Hessian and whole-facet proofs; no source code was copied.

