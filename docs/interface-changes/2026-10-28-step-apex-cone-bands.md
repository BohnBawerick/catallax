## 2026-10-28: band seams for cone faces that end at their apex

Pieces: 7, STEP consumers

Change: A cone face bounded by one vertex-free circular ring and a vertex loop
at its own apex now joins the complete-band encoding of
[the periodic band seams](2026-10-26-periodic-analytic-band-step-seams.md).
Its one `EDGE_LOOP` is the ring, then a longitudinal seam from the apex to the
ring vertex `.F.`, then the same seam `.T.`. The seam's `SEAM_CURVE` carries
pcurves at zero and the rounded full period from the apex's axial height zero,
and the ring keeps its iso-height pcurve. No `VERTEX_LOOP` is written for such a
face. The apex vertex is the source model's vertex. Connected bands share one
meridian as before, so a blind hole's wall and its drill tip are cut together.

This applies to pointed cone primitives, drill-point tips and any other face of
that shape. A file with other faces the band plan does not cover still keeps its
whole legacy encoding. `added_edges` counts the seam. Kernel topology, geometry,
report fields and manifest formats do not change.

Why: The legacy encoding writes the ring with no pcurve and the apex as a vertex
loop. The pinned raw translator then finds no pcurves on the cone or the wall
next to it and makes an apex edge it cannot check, so every blind drill-point
answer failed raw import. This closes the singular-apex gap that the earlier
fragment left open.

Migration: Accept changed STEP bytes for files whose cylinder and cone faces
are all complete bands or apex-ended cone bands and that contain at least one
apex-ended cone band. Update `added_edges` counts for those files. Other files
keep their bytes.
