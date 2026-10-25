## 2026-10-24: retain edge construction names through trim changes

Pieces: topology, boolean, features, application naming and verification tooling

Change: `Model` retains canonical construction origin pairs through boolean edge splits and maximal edge merges. `edge_name` prefers a retained pair that still belongs to opposite incident faces, rather than taking the first origins of already merged faces. `edge_constructions` and `record_edge_constructions` expose this history to topology producers. Exact copies and `map_face_origins` map it with the face origins. The scoreboard adds the P0/seam height edit alongside the existing complete P1/P2/P4 checks.

Why: merging coplanar caps had changed cylinder `b`'s rim name to `a:end` / `b:wall`. Raising cylinder `a` removes that incidence even though cylinder `b`'s construction rim remains.

Migration: topology producers preserve construction pairs when splitting, merging or copying named edges. Pair serialization and reference resolution are unchanged; a pair still needs opposite incident faces. Geometry, export, work allowances and modelling tolerances are unchanged. Timing and edit-survival checks remain development tooling.

