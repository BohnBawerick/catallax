## 2026-10-26: triangulate closed cone and torus bands in layers

Pieces: display mesh, application benchmark

Change: eligible cone and bounded-torus faces with two horizontal sampled rails
join angularly ordered interior layers instead of building periodic CDT copies.
The existing longitude-gap, Hessian and normal-rotation margins select the grid.
Shared vertices remain fixed. Every triangle passes strict UV orientation,
whole-facet chord and rounded-coordinate quality certificates, with exact
fallback for inconclusive interval filters. A failed proposal rolls back all
positions, triangles, face records and unavailable controls before the general
triangulator runs. Nonconstant trims, apices and other periodic domains retain
that triangulator.

Migration: no library API, requested accuracy, resource allowance, model, origin
or STEP changes. Display connectivity can change. The display benchmark now
records `requested_wait_ms` from the start of construction through the completed
mesh, including drawing-bound selection when requested. `display_ms` remains
mesh-only; output serialization remains outside both intervals. The fine-mesh
regression also covers the complete public boss with all four block fillets.
