## 2026-10-25: seed fine revolution meshes from certified shared rails

Pieces: display mesh

Change: eligible cone and bounded-torus interiors use the measured maximum longitude gap in their shared boundary samples. Only the remaining parameter subdivides when that fixed longitude already has a strict chord and angular margin. Sparse trims retain the general grid. Interior columns merge chart-rounding duplicates without modifying shared samples; apex rows retain the existing canonical apex handling. Every returned facet still passes the complete chord and normal certificates.

Why: a balanced cone grid ignores the much finer longitude already present on its boundary, and a bounded torus's nominal grid spacing can overestimate the actual shared-rail gap. Both generate avoidable interior samples and triangulation work.

2026-11-17 trim correction: the general analytic triangulator scales Delaunay distances by the existing grid's U/V step ratio. Mixing angular turns with axial millimetres previously selected long angular diagonals beside sloping trims; centroid refinement then multiplied thin facets. Periodic classification and every chord, angular, edge-length and rounded-facet certificate retain original UVs. Shared boundary samples are unchanged. This is a triangulation metric, not a new geometry approximation or certificate.

A frozen countersink bracket reconstructed from the native fine-mesh fixture, with four 16 mm / 90-degree heads, gives the following release measurements on the same work PC at 0.01 mm chord, no angular target, unchanged 100M allowances and a 10-second mesh deadline. Each timing has one warmup and three measured iterations; these are native mesh timings, not browser edit latency.

| Measure | Before | After |
| --- | ---: | ---: |
| Each side-trimmed countersink, triangles | 58,140 | 5,242 |
| Each untrimmed countersink, triangles | 4,876 | 4,876 |
| Whole bracket, triangles | 132,972 | 27,176 |
| Median certified mesh time, ms | 2,611.45 | 410.93 |

The trimmed-cone aspect bound remains truthful (8,531.71 before, 818.06 after); fixed trim samples can still leave thin facets. No aspect target or tolerance was relaxed. The whole reconstructed fixture is 256 triangles larger than the original Mac report, while reproducing both reported cone counts exactly.

Migration: none. Requested tolerances, user resource limits, exact definitions, STEP output and shared-edge ownership are unchanged. Deterministic display connectivity and vertex counts can change. Rotated-cone regressions compare compact output counts with the measured d0cea699 baseline at identical controls, rather than replacing the temporary-chart budget with an arbitrary output limit.

Prior art: [OCCT's cone range splitter](https://github.com/Open-Cascade-SAS/OCCT/blob/master/src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_ConeRangeSplitter.cxx) scales parameter steps using cone radius and excludes boundary rows. [Shewchuk's anisotropic interpolation discussion](https://people.eecs.berkeley.edu/~jrs/papers/elemj.pdf) separates interpolation accuracy from isotropic element shape. This change keeps Catallax's own Hessian and whole-facet proofs; no source code was copied.

The trim correction follows [OCCT's separation of UV node placement/scaling from final surface-deflection checks](https://dev.opencascade.org/doc/overview/html/occt_user_guides__mesh.html), and [SMLib's separate surface-cell subdivision and trim-polygon clipping](https://github.com/NVIDIA-Omniverse/usd-brep/blob/c9c5979f04cea8c38213004d4f04df0dcfb40d05/source/SMLib/src/NL_Tessellate.cpp#L8675) (`ST_TessSubdivideSrf` / `ST_TessBBoxXTrim`). It reuses Catallax's seed sizes and proofs; no donor code or sampled deviation certificate is introduced.

