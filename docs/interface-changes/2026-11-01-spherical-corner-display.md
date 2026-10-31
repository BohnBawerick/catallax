## 2026-11-01: spherical fillet corners display at an angular target

Pieces: mesh, app

Change: `mesh_body` and `preview_mesh_body` now answer with an angular target for a
sphere face bounded by one loop of three finite great-circle arcs, which is the shape of
a fully filleted corner. Both skip grid seeds inside a boundary segment's diametral disk,
and preview samples those arcs at the certified sampler's chart resolution. Every chord,
rounding and normal proof, the certified edge sampler, the geometry and the public types
are unchanged.

Why: every full spherical fillet corner built, then refused in the app's 15-degree
display with `mesh.precision-or-budget`: a thin triangle on a fixed boundary segment
missed the surface normal by about 30 degrees, and centroid refinement cannot fix that.

Migration: callers that matched that refusal for this family now receive a mesh. Output
that already answered with an angular target inside this family may have different
interior triangles, and preview may sample these arcs and their neighbouring strips more
finely, so compacted indices can differ. Calls with no angular target and every other
face keep their exact output. Invalid options, vertex limits and the spline angular
coverage refusal are unchanged. A preview still makes no accuracy, watertightness or
export claim.
