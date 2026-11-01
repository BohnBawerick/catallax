## 2026-11-01: spherical fillet corners display at the drawing target

Pieces: mesh, app

Change: `mesh_body` and `preview_mesh_body` now answer at the drawing view's 2-degree
target for a sphere face proved to be the positive cone of its three corner rays: one
loop of three finite minor great-circle arcs with the material on the cone's side. Such a
face is meshed in a chart whose poles lie outside the cone, and it takes no grid seed that
an outward interval filter proves is outside the cone. Every chord, rounding, normal,
aspect and closure proof, every cap and work allowance, the geometry and the public types
are unchanged. Preview remains unproved.

Why: each corner evaluated almost a whole sphere of unused seeds, so the body reached
`max_vertices` in its fourth corner's grid. With that waste removed, one skew corner still
kept a thin facet near a chart pole, and centroid refinement moved it along a latitude row
until the subdivision limit refused.

Migration: callers that matched either refusal for this family now receive a mesh. Output
that already answered with an angular target inside this family may change, including its
circular edge polylines and the cylinder strips beside them, so compacted indices can
differ. Calls with no angular target, every sphere face not proved convex, including the
complement or a long-arc trim of the same arcs, and every other face keep their exact
output. Retained handles, exact properties and STEP bytes are unchanged. Invalid options,
vertex limits and the spline angular coverage refusal are unchanged.
